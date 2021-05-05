// Filename main.cpp
// Created on 26 Jul 2016 by Amneet Bhalla
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/CIBMethod.h>
#include <ibamr/CIBMobilitySolver.h>
#include <ibamr/CIBSaddlePointSolver.h>
#include <ibamr/CIBStaggeredStokesSolver.h>
#include <ibamr/DirectMobilitySolver.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/KrylovMobilitySolver.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

//////////////////////////////////////////////////////////////////////////////

struct StructureCtx
{
  std::string name;
  double kappa, eta, rho, vel, mass, R;
  IBTK::Vector X_com, X_com_tether;
  IBTK::Vector U, U_tether;
  IBTK::Vector F;

};// StructureCtx


// Center of mass velocity
void
ConstrainedCOMVel(double /*data_time*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    U_com[1] = 1.0;

    return;
} // ConstrainedCOMOuterVel

void
NetExternalForceTorque(double /*data_time*/, Eigen::Vector3d& F_ext, Eigen::Vector3d& T_ext, void* ctx)
{

    StructureCtx& struct_ctx = *static_cast<StructureCtx*>(ctx);

    F_ext << struct_ctx.F(0), struct_ctx.F(1), 0.0;
    T_ext << 0.0, 0.0, 0.0;

    return;
} // NetExternalForceTorque

void
ConstrainedNodalVel(Vec /*U_k*/, const RigidDOFVector& /*U*/, const Eigen::Vector3d& /*X_com*/, void* /*ctx*/)
{
    // intentionally left blank
    return;
} // ConstrainedNodalVel

ofstream drag_stream;

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 LDataManager* l_data_manager,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
    SAMRAIManager::setMaxNumberPatchDataEntries(2054);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Read default Petsc options
        if (input_db->keyExists("petsc_options_file"))
        {
            std::string petsc_options_file = input_db->getString("petsc_options_file");
            PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, petsc_options_file.c_str(), PETSC_TRUE);
        }

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.

        // INS integrator
        Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        // CIB method
        const unsigned int num_structures = input_db->getIntegerWithDefault("num_structures", 1);
        Pointer<CIBMethod> ib_method_ops =
            new CIBMethod("CIBMethod", app_initializer->getComponentDatabase("CIBMethod"), num_structures);

        // Krylov solver for INS integrator that solves for [u,p,U,L]
        Pointer<CIBStaggeredStokesSolver> CIBSolver =
            new CIBStaggeredStokesSolver("CIBStaggeredStokesSolver",
                                         input_db->getDatabase("CIBStaggeredStokesSolver"),
                                         navier_stokes_integrator,
                                         ib_method_ops,
                                         "SP_");

        // Register the Krylov solver with INS integrator
        navier_stokes_integrator->setStokesSolver(CIBSolver);

        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);

        // Specify structure kinematics
        FreeRigidDOFVector struct_0_free_dofs, struct_1_free_dofs;
        struct_0_free_dofs << 1, 1, 1;
        struct_1_free_dofs << 1, 1, 1;
        ib_method_ops->setSolveRigidBodyVelocity(0, struct_0_free_dofs);
        ib_method_ops->setSolveRigidBodyVelocity(1, struct_1_free_dofs);

        const double K = input_db->getDouble("KAPPA");
	const double eta = input_db->getDouble("ETA");
        const double rho= input_db->getDouble("RHO");
	const double frequency_of_oscillation = input_db->getDouble("FREQUENCY_OF_OSCILLATION");
	const double starting_spring_length = input_db->getDouble("STARTING_SPRING_LENGTH");
        const double max_amp = input_db->getDouble("MAXAMP");
	const double amplitude_of_spring_percent = input_db->getDouble("AMPLITUDE_OF_SPRING_PERCENT");
	const double radius_disc_0 = input_db->getDouble("RADIUS_DISC_0");
	const double radius_disc_1 = input_db->getDouble("RADIUS_DISC_1");

        StructureCtx struct0;
        struct0.name = "disc0";
        struct0.kappa = K;
	struct0.eta = eta;
        struct0.R  = radius_disc_0;
	struct0.U(1) = 0.0;
	struct0.rho = rho;
	struct0.mass = M_PI*struct0.R*struct0.R*struct0.rho;

        StructureCtx struct1;
        struct1.name = "disc1";
        struct1.kappa = K;
	struct1.eta = eta;
        struct1.R  = radius_disc_1;
	struct1.rho = rho;
	struct1.mass = M_PI*struct1.R*struct1.R*struct1.rho;

        ib_method_ops->registerExternalForceTorqueFunction(&NetExternalForceTorque, &struct0, 0);
        ib_method_ops->registerExternalForceTorqueFunction(&NetExternalForceTorque, &struct1, 1);

        ib_method_ops->registerConstrainedVelocityFunction(NULL, &ConstrainedCOMVel, &struct0, 0);
        ib_method_ops->registerConstrainedVelocityFunction(NULL, &ConstrainedCOMVel, &struct1, 1);

        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
            ib_method_ops->registerVisItDataWriter(visit_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Create boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                Pointer<Database> bc_coefs_db = app_initializer->getComponentDatabase(bc_coefs_db_name);
                u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, bc_coefs_db, grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Set physical boundary operator used in spreading.
        ib_method_ops->setVelocityPhysBdryOp(time_integrator->getVelocityPhysBdryOp());

        // Register mobility matrices (if needed)
        std::string mobility_solver_type = input_db->getString("MOBILITY_SOLVER_TYPE");
        if (mobility_solver_type == "DIRECT")
        {
	  int mpi_ranks = SAMRAI_MPI::getNodes();
	  //Need multiple mobility matrices, one for each structure
	  for(int k = 0;k < num_structures; k++)
	    {
	      std::ostringstream mat_name;
	      mat_name << "sphere_mobility_" << k;
	      
	      std::vector<std::vector<unsigned> > struct_ids(1);
	      std::vector<unsigned> prototype_structs(1);

	      //Dens Matrix Type
	      prototype_structs[0] = k; //Type of struct (geo, # IB points etc. Placeholder)
	      struct_ids[0] = prototype_structs; //Actual struct ids of structs

	      //Register the dense matrix with direct solver
	      DirectMobilitySolver* direct_solvers = NULL;
	      CIBSolver->getSaddlePointSolver()->getCIBMobilitySolver()->getMobilitySolvers(NULL, &direct_solvers, NULL);
	      
	      int proc_handling_mob_mat = k % mpi_ranks;
	      direct_solvers->registerMobilityMat(
						  mat_name.str(), prototype_structs,
						  EMPIRICAL, std::make_pair(LAPACK_SVD, LAPACK_SVD), proc_handling_mob_mat);
	      
	      direct_solvers->registerStructIDsWithMobilityMat(mat_name.str(), struct_ids);
	    }
	}

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();

        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
        }
        if (dump_postproc_data)
        {
            output_data(patch_hierarchy,
                        ib_method_ops->getLDataManager(),
                        iteration_num,
                        loop_time,
                        postproc_data_dump_dirname);
        }

        if (SAMRAI_MPI::getRank() == 0)
        {
            drag_stream.open("./Lambda/drag.txt", std::ios_base::out | ios_base::trunc);
            drag_stream.precision(10);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
	double pastU, pastPosX,pastPosY;
	double current_spring_length = 0.0;
	double dist_bw_discs[2];
	const double amplitude_of_spring = amplitude_of_spring_percent*max_amp;
	double desired_spring_length = starting_spring_length;
	double direction_vector[2];

	//ofstream force_data;
	//force_data.open("fd.txt",ios::out | ios::trunc);
	ofstream position_data;
	position_data.open("pd.txt", ios::out | ios::trunc);

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();

            pout << "Advancing hierarchy by timestep size dt = " << dt << "\n";
            if (time_integrator->atRegridPoint()) navier_stokes_integrator->setStokesSolverNeedsInit();
            if (ib_method_ops->flagRegrid())
            {
                time_integrator->regridHierarchy();
                navier_stokes_integrator->setStokesSolverNeedsInit();
            }

            const Eigen::Vector3d& X_com_0 = ib_method_ops->getCurrentBodyCenterOfMass(0);
            const Eigen::Vector3d& X_com_1 = ib_method_ops->getCurrentBodyCenterOfMass(1);

            RDV U0, U1;
            ib_method_ops->getCurrentRigidBodyVelocity(0, U0);
            ib_method_ops->getCurrentRigidBodyVelocity(1, U1);

	    pout << "RigidBodyVelocity1 = " << U1(1) << "\n";

	    //Store Past position for active disk
	    pastPosX = struct1.X_com(0);
	    pastPosY = struct1.X_com(1);

	    //Store Position and Velocity Data of Discs
            for (int d = 0; d < NDIM; ++d)
            {
                struct0.X_com(d) = X_com_0(d);
		struct0.X_com_tether(d) = X_com_1(d);
                struct0.U(d) = U0(d);
		struct0.U_tether(d) = U1(d);

                struct1.X_com(d) = X_com_1(d);
		struct1.X_com_tether(d) = X_com_0(d);
                struct1.U(d) = U1(d);
		struct1.U_tether(d) = U0(d);

		dist_bw_discs[d] = X_com_0(d) - X_com_1(d); //dist (in each dim) bw discs
            }

	    //Calculate current/desired spring lengths
	    current_spring_length = sqrt(dist_bw_discs[0]*dist_bw_discs[0] + dist_bw_discs[1]*dist_bw_discs[1]);
	    desired_spring_length = starting_spring_length + amplitude_of_spring*(sin(2.0*M_PI*frequency_of_oscillation*(loop_time+dt)));

	    //Unit Vector of spherobot system
	    direction_vector[0] = dist_bw_discs[0]/current_spring_length;
	    direction_vector[1] = dist_bw_discs[1]/current_spring_length;

	    //Step1: Calculate Forces applied to each disc
	    struct0.F(0) = struct0.kappa*(0.0 - X_com_0(0)); //constrains x-pos to zero
	    struct0.F(1) = struct0.kappa*0.5*(desired_spring_length - current_spring_length);

	    struct1.F(0) = struct1.kappa*(0.0 - X_com_1(0));
	    struct1.F(1) = -1.0*struct1.kappa*0.5*(desired_spring_length - current_spring_length);

	    /*//Calculate current velocity of active disk
	    double dist = ((struct1.X_com(0) - pastPosX)*(struct1.X_com(0) - pastPosX) +
			   (struct1.X_com(1) - pastPosY)*(struct1.X_com(1) - pastPosY));
	    pastU = sqrt(dist)/dt;
	    if(loop_time == 0.0)
	      pout << "Current Velocity1 = " << pastU << "\n";
	    
	    //Step1
	    if(struct1.U(0) <= struct1.vel*cos(theta)){
	      struct1.F(0) = struct1.kappa*struct1.vel*dt*cos(theta); //Active sphere. Moves at an angle theta
	    }
	    else{
	      struct1.F(0) = 0.0;
	    }
	    if(struct1.U(1) <= struct1.vel*sin(theta)){
	      struct1.F(1) = struct1.kappa*struct1.vel*dt*sin(theta);
	    }
	    else{
	      struct1.F(1) = 0.0;
	    }
	    
	    struct0.F(0) = 0.0; //Passive Sphere
	    struct0.F(1) = 0.0;*/
	    
	    pout << "Done Finding Forces! Advance Time Hierarchy!\n";
	    
	    //Steps 2,3,4
            time_integrator->advanceHierarchy(dt);
	    pout << "Hierarchy Advanced!\n";
            loop_time += dt;

	    //OUTPUT
	    position_data << X_com_0(0) << " " << X_com_0(1) << " " << X_com_1(0) << " " << X_com_1(1) << " " << current_spring_length << " " << desired_spring_length << " " << loop_time <<"\n" << std::flush;
	    //force_data << struct0.F(0) << " " << struct0.F(1) << " " << struct1.F(0) << " " << struct1.F(1) << " " << loop_time << "\n" << std::flush;

	    pout << "The position of circle 0 is " << X_com_0(0) << " " << X_com_0(1) << "\n";
	    pout << "The position of circle 1 is " << X_com_1(0) << " " << X_com_1(1) << "\n";
            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                output_data(patch_hierarchy,
                            ib_method_ops->getLDataManager(),
                            iteration_num,
                            loop_time,
                            postproc_data_dump_dirname);
            }
        }
	
	//force_data.close();
	position_data.close();

         if (SAMRAI_MPI::getRank() == 0)
         {
             drag_stream.close();
         }
          
        // Cleanup boundary condition specification objects (when necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main

void
output_data(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
            LDataManager* /*l_data_manager*/,
            const int iteration_num,
            const double loop_time,
            const string& /*data_dump_dirname*/)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    return;
} // output_data
