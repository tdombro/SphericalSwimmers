#include "update_springs.h"
#include <ibamr/IBSpringForceSpec.h>

void
update_springs(
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const l_data_manager,
    const double current_time,
    const double dt)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();

    static const double pi = 4*atan(1);
    static const double L1 = 0.025; // length of computational domain (meters)
    static const int Nbots = 100; //Number of swimmers
    static const int N1 = 256; // number of cartesian grid meshwidths at the finest level of the AMR grid
    static const int nvert = 26; //number of vertices in skeleton
    static const double dX = L1/(1.0*N1);
    static const double amp = 0.002*0.8;
    static const double upper_radius = 0.002 - 1.25*dX;
    static const double lower_radius = 0.001 - 1.25*dX;
    static const double upper_half_radius = upper_radius/2.0;
    static const double lower_half_radius = lower_radius/2.0;

    double xLength, yLength;

    // Find out the Lagrangian index ranges.
    const std::pair<int,int>& wing_lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);

    // Get the LMesh (which we assume to be associated with the finest level of
    // the patch hierarchy).  Note that we currently need to update both "local"
    // and "ghost" node data.
    Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
    vector<LNode*> nodes;
    nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
    nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());

    // Update the spring lengths in their associated spring specs.
    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);
    for (vector<LNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
      {
        LNode* node_idx = *it;
        IBSpringForceSpec* spring_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
		
        if (spring_spec == NULL) continue;  // skip to next node

        // Here we update the resting length of the spring
        //
        // NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the number of Lagrangian points)
        //        resting_length    	is the resting length of the current Lagrangian point that is the "master index'
	//		  						Since there may be more than one spring associated with each point, it's a vector.
        //        resting_length[0]  	is the resting length of the first spring
        //        resting_length[1]  	would be the resting length of the second spring associated with that Lagrangian point.
        //
        // In this example, the resting length is increased by 0.01*dt each time step.

        const int lag_idx = node_idx->getLagrangianIndex();
	//std::vector<double>& resting_length = spring_spec->getRestingLengths();
	//The above is the old way of doing this.
	//double resting_length = spring_spec->getParameters()[0][1];
	std::vector<double>& spring0 = spring_spec->getParameters()[0];
	
	//Note that you can also getStiffnesses
	//double spring_stiffness = spring_spec->getParameters()[0][0];
	//resting_length=0.75 + 0.5*0.6*0.8*sin(2.0*M_PI*8.0*(current_time+dt));
	for(int i=0;i<Nbots;i++){
	  if (lag_idx == 6 + nvert*i){
	    pout << "b4: lag_idx = " << lag_idx << "\tresting_length = " << spring0[1] << "\n";
	    pout << "resting_length contains " << spring0.size() << " elements \n";
	    spring0[1] = 0.005 + amp*sin(2.0*M_PI*10.0*(current_time+dt));
	    pout << "a4: resting_length = " << spring0[1] << "\n";
	  }
	  
	  if(lag_idx == 2 + nvert*i){
	    pout << "b4: lag_idx = " << lag_idx << "\tresting_length = " << spring0[1] << "\n";
	    pout << "resting_length contains " << spring0.size() << " elements \n";
	    xLength = upper_half_radius - lower_radius;
	    yLength = 0.005 + amp*sin(2.0*M_PI*10.0*(current_time+dt));
	    spring0[1] = sqrt(xLength*xLength + yLength*yLength);
	    pout << "a4: resting_length = " << spring0[1] << "\n";
	    std::vector<double>& spring1 = spring_spec->getParameters()[1];
	    pout << "b4: lag_idx = " << lag_idx << "\tresting_length = " << spring1[1] << "\n";
	    xLength = upper_half_radius + lower_radius;
	    yLength = 0.005 + amp*sin(2.0*M_PI*10.0*(current_time+dt));
	    spring1[1] = sqrt(xLength*xLength + yLength*yLength);
	    pout << "a4: resting_length = " << spring1[1] << "\n";
	  }

	  if(lag_idx == 10 + nvert*i){
	    pout << "b4: lag_idx = " << lag_idx << "\tresting_length = " << spring0[1] << "\n";
	    pout << "resting_length contains " << spring0.size() << " elements \n";
	    xLength = upper_half_radius - lower_radius;
	    yLength = 0.005 + amp*sin(2.0*M_PI*10.0*(current_time+dt));
	    spring0[1] = sqrt(xLength*xLength + yLength*yLength);
	    pout << "a4: resting_length = " << spring0[1] << "\n";
	    std::vector<double>& spring1 = spring_spec->getParameters()[1];
	    pout << "b4: lag_idx = " << lag_idx << "\tresting_length = " << spring1[1] << "\n";
	    xLength = upper_half_radius + lower_radius;
	    yLength = 0.005 + amp*sin(2.0*M_PI*10.0*(current_time+dt));
	    spring1[1] = sqrt(xLength*xLength + yLength*yLength);
	    pout << "a4: resting_length = " << spring1[1] << "\n";
	  }

	}

	/*if(lag_idx == 2){
	  pout << "b4: lag_idx = " << lag_idx << "\tresting_length = " << spring0[1] << "\n";
	  pout << "resting_length contains " << spring0.size() << " elements \n";
	  spring0[1] = 0.375 + 0.5*0.3*0.8*sin(2.0*M_PI*8.0*(current_time+dt));
	  pout << "a4: resting_length = " << spring0[1] << "\n";
	}
	
	if(lag_idx == 13){
	  pout << "b4: lag_idx = " << lag_idx << "\tresting_length = " << spring0[1] << "\n";
	  pout << "resting_length contains " << spring0.size() << " elements \n";
	  xLength = upper_half_radius - lower_radius;
	  yLength = 0.375 + 0.5*0.3*0.8*sin(2.0*M_PI*8.0*(current_time+dt));
	  spring0[1] = sqrt(xLength*xLength + yLength*yLength);
	  pout << "a4: resting_length = " << spring0[1] << "\n";
	  std::vector<double>& spring1 = spring_spec->getParameters()[1];
	  pout << "b4: lag_idx = " << lag_idx << "\tresting_length = " << spring1[1] << "\n";
	  xLength = upper_half_radius + lower_radius;
	  yLength = 0.375 + 0.5*0.3*0.8*sin(2.0*M_PI*8.0*(current_time+dt));
	  spring1[1] = sqrt(xLength*xLength + yLength*yLength);
	  pout << "a4: resting_length = " << spring1[1] << "\n";
	}

	if(lag_idx == 14){
	  pout << "b4: lag_idx = " << lag_idx << "\tresting_length = " << spring0[1] << "\n";
	  pout << "resting_length contains " << spring0.size() << " elements \n";
	  xLength = upper_half_radius - lower_radius;
	  yLength = 0.375 + 0.5*0.3*0.8*sin(2.0*M_PI*8.0*(current_time+dt));
	  spring0[1] = sqrt(xLength*xLength + yLength*yLength);
	  pout << "a4: resting_length = " << spring0[1] << "\n";
	  std::vector<double>& spring1 = spring_spec->getParameters()[1];
	  pout << "b4: lag_idx = " << lag_idx << "\tresting_length = " << spring1[1] << "\n";
	  xLength = upper_half_radius + lower_radius;
	  yLength = 0.375 + 0.5*0.3*0.8*sin(2.0*M_PI*8.0*(current_time+dt));
	  spring1[1] = sqrt(xLength*xLength + yLength*yLength);
	  pout << "a4: resting_length = " << spring1[1] << "\n";
	  }*/

      }
    return;
}// update_springs
