// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////////////// INCLUDES ///////////////////////////////////////////

#include "ibamr/namespaces.h"

#include "Box.h"
#include "ForceProjector.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"
#include "tbox/PIO.h"


namespace IBTK
{
void
callForceProjectorCallBackFunction(const double current_time, const double new_time, const int cycle_num, void* ctx)
{
    static ForceProjector* ptr_forceprojector = static_cast<ForceProjector*>(ctx);
    if (cycle_num == 0)
    {
        ptr_forceprojector->calculateLagrangianBodyForce(new_time, current_time);
        ptr_forceprojector->calculateEulerianBodyForce(new_time, current_time);
    }

    return;

} // callForceProjectorCallBackFunction

ForceProjector::ForceProjector(const std::string& object_name,
                               LDataManager* lag_data_manager,
                               Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                               Pointer<Database> input_db,
                               const std::string solver_type)
    : d_object_name(object_name),
      d_lag_data_manager(lag_data_manager),
      d_patch_hierarchy(patch_hierarchy),
      d_solver_type(solver_type),
      d_circle_center(NDIM)
{
    d_n_large = d_n_small = 1;
    const IntVector<NDIM> ib_ghosts(3);// = getMinimumGhostCellWidth();
    // Initialize  variables & variable contexts associated with Eulerian forces.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_body_force_context = var_db->getContext(d_object_name + "::BODYFORCE");
    if (d_solver_type == "STAGGERED")
        d_body_force_var = new SideVariable<NDIM, double>(d_object_name + "::BodyForce_sc_var");
    if (d_solver_type == "COLLOCATED")
        d_body_force_var = new CellVariable<NDIM, double>(d_object_name + "::BodyForce_cc_var", NDIM);
    d_body_force_idx = var_db->registerVariableAndContext(d_body_force_var, d_body_force_context, ib_ghosts);

    getFromInput(input_db);

    return;

} // ForceProjector

ForceProjector::~ForceProjector()
{
    // intentionally left blank
    return;

} //~ForceProjector

void
ForceProjector::getFromInput(Pointer<Database> input_db)
{
    d_radius_0 = input_db->getDouble("radius_0");
    d_radius_1 = input_db->getDouble("radius_1");
    d_delta    = input_db->getDouble("interaction_range");
    d_kappa    = input_db->getDouble("kappa");

    d_wall_left_x   = input_db->getDouble("wall_left");
    d_wall_right_x  = input_db->getDouble("wall_right");
    d_wall_up_y     = input_db->getDouble("wall_up");
    d_wall_bottom_y = input_db->getDouble("wall_bottom");
    d_enable_wall   = input_db->getBool("enable_wall");

    d_enable_circle_bdry = input_db->getBoolWithDefault("enable_circle_bdry", false);
    d_circle_bdry_radius = input_db->getDoubleWithDefault("circle_bdry_radius", 0);
    input_db->getDoubleArray("circle_center", &d_circle_center[0], NDIM);

    return;
} // getFromInput

void
ForceProjector::registerLagrangianQuantityName(const std::string& lag_quantity_name)
{
    registerLagrangianQuantitiesName(std::vector<std::string>(1, lag_quantity_name));

    return;

} // registerLagrangianQuantityName

void
ForceProjector::registerLagrangianQuantitiesName(const std::vector<std::string>& lag_quantities_name)
{
    const unsigned size = lag_quantities_name.size();
    for (unsigned i = 0; i < size; ++i)
    {
        d_lag_quantities_name.push_back(lag_quantities_name[i]);
    }

    return;

} // registerLagrangianQuantitiesName

void
ForceProjector::associateVolumeElement(const double vol_lag_pt)
{
    d_vol_lag_pt = vol_lag_pt;

    return;

} // associateVolumeElement

void
ForceProjector::calculateLagrangianBodyForce(const double /*new_time*/, const double current_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    d_lag_force.clear();
    d_lag_force.resize(finest_ln + 1);

    int lag_pts_disc = 0;

    //#ifdef MYLOG
    //printf("ForceProjector::calculateLagrangianBodyForce::\n");
    //#endif

    for (unsigned int sid = 0; sid < d_struct_pos_vec.size(); ++sid){
        std::vector<double> conf = d_struct_pos_vec[sid];
        //#ifdef MYLOG
        //printf("\tstruct %d [%g, %g, %g, %g, %g, %g]\n", sid, conf[0], conf[1], conf[2], conf[3], conf[4], conf[5]);
        //#endif
        lag_pts_disc += conf[4] + conf[5];
    }

    //#ifdef MYLOG
    //printf("\ttotal lag pts disc = %d\n", lag_pts_disc);
    //#endif

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_lag_data_manager->levelContainsLagrangianData(ln)) continue;
        d_lag_force[ln] = d_lag_data_manager->createLData(d_object_name + "::lag_force_data", ln, NDIM, false);
    }

    /*Pointer<LData> ptr_x_lag_data_current(nullptr);
    ptr_x_lag_data_current = d_lag_data_manager->getLData("X", finest_ln);
    const boost::multi_array_ref<double, 2>& X_data_current = *ptr_x_lag_data_current->getLocalFormVecArray();
    */

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_lag_data_manager->levelContainsLagrangianData(ln)) continue;


        std::vector<int> struct_ids = d_lag_data_manager->getLagrangianStructureIDs(ln);
        for (unsigned int sid = 0; sid < struct_ids.size(); ++sid){
            //std::pair<int,int> 	idx_range = d_lag_data_manager->getLagrangianStructureIndexRange(struct_ids[sid], ln);

            //#ifdef MYLOG
            //pout << "ForceProjector::calculateLagrangianBodyForce:: struct " << struct_ids[sid] << ", lag idx range = " << idx_range.first << ", " << idx_range.second << "\n";
            //#endif
        }

        // Get ponter to LData corresponding to lagrangian force.
        boost::multi_array_ref<double, 2>& F_data = *d_lag_force[ln]->getLocalFormVecArray();
        const Pointer<LMesh> mesh = d_lag_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        //int total_nodes = 0;

        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int local_idx = node_idx->getLocalPETScIndex();
            double* const F = &F_data[local_idx][0];
            int lag_idx = node_idx->getLagrangianIndex();


            //const double* const X_current = &X_data_current[local_idx][0];
            //total_nodes++;
            /*if(lag_idx == 0){
              printf("ForceProjector::calculateLagrangianBodyForce HONG: called at  t = %g, X[%d] = [%g, %g]\n",
                current_time, local_idx, X_current[0], X_current[1] );
            }
            */
            for (int d = 0; d < 2; ++d){
                if (lag_idx < lag_pts_disc)
                    F[d] = d_repul_force[lag_idx][d]; // * d_vol_lag_pt;
                else
                    F[d] = 0.0;
            }

            //#ifdef MYLOG
            /*if (lag_idx == d_n_large + d_n_small - 1)
                printf("ForceProjector::calculateLagrangianBodyForce: time %g, F_small = [%g, %g], lag_idx %d\n",
                  current_time, F[0], F[1], lag_idx);

            if (lag_idx == d_n_large - 1)
                printf("ForceProjector::calculateLagrangianBodyForce: time %g, F_large = [%g, %g], lag_idx %d\n",
                        current_time, F[0], F[1], lag_idx);
            #endif*/
        }

        d_lag_force[ln]->restoreArrays();
    } // all levels

    return;
} // calculateLagrangianBodyForce

void
ForceProjector::calculateEulerianBodyForce(const double /*new_time*/, const double current_time)
{
    // allocate patch data for Eulerian forcing.
    const int coarsest_ln = 0;
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_body_force_idx)) level->deallocatePatchData(d_body_force_idx);
        level->allocatePatchData(d_body_force_idx, current_time);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            if (d_solver_type == "STAGGERED")
            {
                Pointer<SideData<NDIM, double> > body_force_data = patch->getPatchData(d_body_force_idx);
                body_force_data->fill(0.0);
            }
            else if (d_solver_type == "COLLOCATED")
            {
                Pointer<CellData<NDIM, double> > body_force_data = patch->getPatchData(d_body_force_idx);
                body_force_data->fill(0.0);
            }
            else
            {
                TBOX_ERROR("ForceProjector::calculateEulerianBodyForce() "
                           << "UNKNOWN SOLVER ENCOUNTERED" << std::endl);
            }

        } // iterate over patches

    } // all levels.

    // spread the lagrangian force from finest level to the finest level.
    std::vector<Pointer<LData> > F_data(finest_ln + 1, Pointer<LData>(NULL));
    std::vector<Pointer<LData> > X_data(finest_ln + 1, Pointer<LData>(NULL));

    // Fill in the above vectors at the finest level.
    F_data[finest_ln] = d_lag_force[finest_ln];
    X_data[finest_ln] = d_lag_data_manager->getLData("X", finest_ln);

    // Spread the deformation velocities.
    d_lag_data_manager->spread(d_body_force_idx, F_data, X_data, (RobinPhysBdryPatchStrategy*)NULL);

    return;

} // calculateEulerianBodyForce



void
ForceProjector::updateStructCoordinates(const std::vector<std::vector<double>> & struct_pos_vec)
{
    d_struct_pos_vec.clear();

    for (unsigned int sid = 0; sid < struct_pos_vec.size(); ++ sid){
        d_struct_pos_vec.push_back(struct_pos_vec[sid]);
    }

    calcRepulsiveForces();
    return;
} //

void
ForceProjector::calcRepulsiveForces()
{
    int num_struct = d_struct_pos_vec.size();
    Eigen::Vector2d ll, ls, sl, ss, f_small[num_struct], f_large[num_struct];

    // reset force
    for (int i = 0; i < num_struct; ++i) {
        f_small[i].setZero();
        f_large[i].setZero();
    }

    // repulsive force btw bots
    for (int i = 0; i < num_struct - 1; ++i){
        Eigen::Vector2d largeI (d_struct_pos_vec[i][0], d_struct_pos_vec[i][1]);
        Eigen::Vector2d smallI (d_struct_pos_vec[i][2], d_struct_pos_vec[i][3]);

        for (int j = i + 1; j < num_struct; ++j){
            Eigen::Vector2d largeJ (d_struct_pos_vec[j][0], d_struct_pos_vec[j][1]);
            Eigen::Vector2d smallJ (d_struct_pos_vec[j][2], d_struct_pos_vec[j][3]);

            ll = largeI - largeJ;
            ls = largeI - smallJ;

            sl = smallI - largeJ;
            ss = smallI - smallJ;

            // large - large
            if ( ll.norm() < (2*d_radius_0 + d_delta ) ){
                double r = ll.norm() - (2*d_radius_0 + d_delta);
                f_large[i] +=  d_kappa * ll * r * r;
                f_large[j] += -d_kappa * ll * r * r;
            }

            // large small
            if ( ls.norm() < (d_radius_0 + d_radius_1 + d_delta ) ){
                double r = ls.norm() - (d_radius_0 + d_radius_1 + d_delta);
                f_large[i] +=  d_kappa * ls * r * r;
                f_small[j] += -d_kappa * ls * r * r;
            }

            // small large
            if ( sl.norm() < (d_radius_0 + d_radius_1 + d_delta ) ){
                double r = sl.norm() - (d_radius_0 + d_radius_1 + d_delta);
                f_small[i] +=  d_kappa * sl * r * r;
                f_large[j] += -d_kappa * sl * r * r;
            }

            // small small
            if ( ss.norm() < (d_radius_1 + d_radius_1 + d_delta ) ){
                double r = ss.norm() - (d_radius_1 + d_radius_1 + d_delta);
                f_small[i] +=  d_kappa * ss * r * r;
                f_small[j] += -d_kappa * ss * r * r;
            }
        }// for j
    } // for i

    // repulsive from wall
    // functional form of this force is the same as in bot-bot repulsion
    if (d_enable_wall){
      Eigen::Vector2d vc;
      double r;

      for (int i = 0; i < num_struct; ++i){
          Eigen::Vector2d large (d_struct_pos_vec[i][0], d_struct_pos_vec[i][1]);
          Eigen::Vector2d small (d_struct_pos_vec[i][2], d_struct_pos_vec[i][3]);

          // left_wall
          Eigen::Vector2d large_left (2*d_wall_left_x - large(0), large(1));
          Eigen::Vector2d small_left (2*d_wall_left_x - small(0), small(1));

          vc = large - large_left;
          r = vc.norm() - 2*d_radius_0 - d_delta;
          if ( r < 0 ) f_large[i] += d_kappa * vc * r * r;

          vc = small - small_left;
          r = vc.norm() - 2*d_radius_1 - d_delta;
          if (r < 0) f_small[i] += d_kappa * vc * r * r;

          // right_wall
          Eigen::Vector2d large_right (2*d_wall_right_x - large(0), large(1));
          Eigen::Vector2d small_right (2*d_wall_right_x - small(0), small(1));
          vc = large - large_right;
          r = vc.norm() - 2*d_radius_0 - d_delta;
          if (r < 0) f_large[i] += d_kappa * vc * r * r;


          vc = small - small_right;
          r = vc.norm() - 2*d_radius_1 - d_delta;
          if (r < 0) f_small[i] += d_kappa * vc * r * r;


          // wall_bottom
          Eigen::Vector2d large_bottom (large(0), 2*d_wall_bottom_y - large(1));
          Eigen::Vector2d small_bottom (small(0), 2*d_wall_bottom_y - small(1));
          vc = large - large_bottom;
          r = vc.norm() - 2*d_radius_0 - d_delta;
          if (r < 0) f_large[i] += d_kappa * vc * r * r;

          vc = small - small_bottom;
          r = vc.norm() - 2*d_radius_1 - d_delta;
          if (r < 0) f_small[i] += d_kappa * vc * r * r;

          // wall_up
          Eigen::Vector2d large_up (large(0), 2*d_wall_up_y - large(1));
          Eigen::Vector2d small_up (small(0), 2*d_wall_up_y - small(1));
          vc = large - large_up;
          r = vc.norm() - 2*d_radius_0 - d_delta;
          if (r < 0) f_large[i] += d_kappa * vc * r * r;

          vc = small - small_up;
          r = vc.norm() - 2*d_radius_1 - d_delta;
          if (r < 0) f_small[i] += d_kappa * vc * r * r;

        }// num_struct
    }


    // circle interaction
    if (d_enable_circle_bdry){
      // center of the circle
      Eigen::Vector2d rC(d_circle_center[0], d_circle_center[1]);
      Eigen::Vector2d CD, CS, rS, rI, vc; // S on disc surface, I is image of disc D, C is center of circle
      double r;

      // loop over all struct
      for (int i = 0; i < num_struct; ++i){
          for (int discId = 0; discId < 2; ++discId){
              Eigen::Vector2d rD(d_struct_pos_vec[i][2*discId], d_struct_pos_vec[i][2*discId + 1]);
              CD = rD - rC;
              CS = CD * d_circle_bdry_radius / CD.norm();
              rS = rC + CS;
              rI = 2*rS - rD;
              vc = rD - rI;

              // large
              if (discId == 0){
                  r = vc.norm() - 2*d_radius_0 - d_delta;
                  if (r < 0) f_large[i] += d_kappa * vc * r * r;
              }/*small*/
              else{
                  r = vc.norm() - 2*d_radius_1 - d_delta;
                  if (r < 0) f_small[i] += d_kappa * vc * r * r;
              }
          }//disc ID
      }// struct ID
    }// circle


    // map into lag indices
    d_repul_force.clear();
    int tot_lag_idx = 0;
    for (int sid = 0; sid < num_struct; ++sid){
        std::vector<double> f;

        int n_large_idx = d_struct_pos_vec[sid][4];
        int n_small_idx = d_struct_pos_vec[sid][5];

        d_n_large = n_large_idx;
        d_n_small = n_small_idx;

        tot_lag_idx += n_large_idx + n_small_idx;

        // large disc
        f.clear();
        f.push_back(f_large[sid](0) / n_large_idx);
        f.push_back(f_large[sid](1) / n_large_idx);
        for (int i = 0; i < n_large_idx; ++i) d_repul_force.push_back(f);

        // small disck
        f.clear();
        f.push_back(f_small[sid](0) / n_small_idx);
        f.push_back(f_small[sid](1) / n_small_idx);
        for (int i = 0; i < n_small_idx; ++i) d_repul_force.push_back(f);
    }

    /*#ifdef MYLOG
    printf("ForceProjector::calcRepulsiveForces:: tot lag pts %d, large [0, %d], small = [%d, %d]\n",
      tot_lag_idx, d_n_large -1, d_n_large, d_n_large + d_n_small - 1 );
    #endif*/

    return;
} // calcRepulsiveForces

} // namespace IBTK
