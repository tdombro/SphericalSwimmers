#ifndef included_update_target_point_positions
#define included_update_target_point_positions

#include <ibamr/app_namespaces.h>
#include <ibtk/LDataManager.h>

#include <ibtk/LData.h> //Added for MODULE 2016-11

/*
 * Update the positions of the target point specifications.
 */
void
update_target_point_positions(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    LDataManager* l_data_manager,
    const double current_time,
    const double dt);

#endif //#ifndef included_update_target_point_positions
