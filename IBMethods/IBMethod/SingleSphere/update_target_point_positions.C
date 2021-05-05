#include "update_target_point_positions.h"
#include <ibamr/IBTargetPointForceSpec.h>

void
update_target_point_positions(
			      tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
			      LDataManager* const l_data_manager,
			      const double current_time,
			      const double dt)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();

    static const double pi = 4*atan(1);
    //static const double V = 0.2; //velocity of the wing during translation (meters/sec)
    //static const double time_to_acc = 0.1; //time to accelerate up to V (s)

    static const double L1 = 16.0; // length of computational domain (meters)
    static const int N1 = 400; // number of cartesian grid meshwidths at the finest level of the AMR grid

    static const double KC  = 5.0;
    static const double RADIUS = 0.5;
    static const double DIAMETER = 2.0*RADIUS;
    static const double FREQUENCY = 10.0;
    static const double U_0       = KC*FREQUENCY*DIAMETER;
    static const double OMEGA = 2.0*pi*FREQUENCY;

    // Find out the Lagrangian index ranges.
    const std::pair<int,int>& wing_lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);

    // Get the LMesh (which we assume to be associated with the finest level of
    // the patch hierarchy).  Note that we currently need to update both "local"
    // and "ghost" node data.
    Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
    vector<LNode*> nodes;
    nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
    nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());

    // Update the target point positions in their associated target point force
    // specs.
    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);
    for (vector<LNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            LNode* node_idx = *it;
            IBTargetPointForceSpec* force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
            if (force_spec == NULL) continue;  // skip to next node

            // Here we update the position of the target point.
            //
            // NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the number of Lagrangian points)
            //        X_target     is the target position of the target point
            //        X_target[0]  is the x component of the target position
            //        X_target[1]  is the y component of the target position
            //        X_target[2]  is the z component of the target position (for a 3D simulation)
            //
            // The target position is shifted to the left or right by the
            // increment dt*V

            const int lag_idx = node_idx->getLagrangianIndex();
            //Depending on the version of IBAMR, you need to select one of the ways of accessing target point positions
            //TinyVector<double,NDIM>& X_target = force_spec->getTargetPointPosition();
            //IBTK::Vector<double,NDIM>& X_target = force_spec->getTargetPointPosition();
            Point& X_target = force_spec->getTargetPointPosition();
            //Print CM target pt before
            if(lag_idx == 0)
                pout << "Bot1up: b4: lag_idx = " << lag_idx << "\txpos = " << X_target[0] << "\typos = " << X_target[1] << "\n";
            //Move target pts by V*dt for all of the cylinder
            X_target[1] += U_0*cos(OMEGA*current_time)*dt;
            //print CM target pt after
            if(lag_idx == 0)
                pout << "a4: lag_idx = " << lag_idx << "\txpos = " << X_target[0] << "\typos = " << X_target[1] << "\n";

    }
  return;
}// update_target_point_positions
