// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

//////////////////////////// INCLUDES /////////////////////////////////////////
#include "ibamr/namespaces.h"

#include "CartesianPatchGeometry.h"
#include "Kinematics.h"
#include "PatchLevel.h"
#include "tbox/MathUtilities.h"
#include "tbox/SAMRAI_MPI.h"

#include "muParser.h"

#include <cmath>
#include <fstream>
#include <iostream>

namespace IBAMR
{
namespace
{
static const double PII = 3.1415926535897932384626433832795;
} // namespace

///////////////////////////////////////////////////////////////////////

Kinematics::Kinematics(const std::string& object_name,
                                 Pointer<Database> input_db,
                                 LDataManager* l_data_manager,
                                 Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                 bool register_for_restart)
    : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
      d_current_time(0.0),
      d_kinematics_vel(NDIM),
      d_shape(NDIM),
      d_center_of_mass(3),
      d_incremented_angle_from_reference_axis(3),
      d_tagged_pt_position(3),
      d_mesh_width(NDIM)
{
    d_radius_0 = input_db->getDouble("radius_0");
    d_radius_1 = input_db->getDouble("radius_1");
    d_center_to_center = input_db->getDouble("center_to_center");
    d_amplitude = input_db->getDouble("amplitude");
    d_frequency = input_db->getDouble("frequency");
    d_large_lag_pts = input_db->getInteger("num_lag_pts_large");

    // set how the immersed body is layout in reference frame.
    setImmersedBodyLayout(patch_hierarchy);

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;

} // Kinematics

Kinematics::~Kinematics()
{
    return;

} // ~Kinematics

void
Kinematics::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_current_time", d_current_time);
    db->putDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->putDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->putDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;

} // putToDatabase

void
Kinematics::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }

    d_current_time = db->getDouble("d_current_time");
    db->getDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->getDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->getDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;
} // getFromRestart

void
Kinematics::setImmersedBodyLayout(Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
{
    printf("Kinematics::setImmersedBodyLayout:: Called\n");

    const StructureParameters& struct_param = getStructureParameters();
    const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    const int finest_ln = struct_param.getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln == finest_ln);
    const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();
    const int total_lag_pts = idx_range[0].second - idx_range[0].first;

    for (int d = 0; d < NDIM; ++d)
    {
        d_kinematics_vel[d].resize(total_lag_pts);
        d_shape[d].resize(total_lag_pts);
    }

    // Get Background mesh related data.
    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(finest_ln);
    PatchLevel<NDIM>::Iterator p(level);
    Pointer<Patch<NDIM> > patch = level->getPatch(p());
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    for (int dim = 0; dim < NDIM; ++dim)
    {
        d_mesh_width[dim] = dx[dim];
    }

    // No. of points on the backbone and till head.
    const int Disc_0_Nx = static_cast<int>(floor(d_radius_0 / d_mesh_width[0]));
    const int Disc_1_Nx = static_cast<int>(floor(d_radius_1 / d_mesh_width[0]));

    d_total_lag_pts = total_lag_pts;
    d_ImmersedBodyData.clear();
    for (int i = -Disc_0_Nx; i <= Disc_0_Nx; ++i)
    {
        const double x = i * d_mesh_width[0];
        const double section = sqrt(d_radius_0 * d_radius_0 - x * x);
        const int NumPtsInSection = static_cast<int>(floor(section / d_mesh_width[1]));
        d_ImmersedBodyData.insert(std::make_pair(x, NumPtsInSection));
    }

    for (int i = -Disc_1_Nx; i <= Disc_1_Nx; ++i)
    {
        const double x =  i * d_mesh_width[0];
        const double section = sqrt(d_radius_1 * d_radius_1 - x * x);;
        const int NumPtsInSection = static_cast<int>(floor(section / d_mesh_width[1]));
        d_ImmersedBodyData.insert(std::make_pair(x + d_center_to_center, NumPtsInSection));
    }

    return;

} // setImmersedBodyLayout

void
Kinematics::setSpecificVelocity(const double time,
                                        const std::vector<double>& incremented_angle_from_reference_axis,
                                        const std::vector<double>& center_of_mass,
                                        const std::vector<double>& tagged_pt_position)
{
    double t  = time;
    double vc = d_amplitude * 2 * PII * d_frequency * std::cos(2* PII* d_frequency * t);
    double rc3 = d_radius_0*d_radius_0*d_radius_0 + d_radius_1*d_radius_1*d_radius_1;
    //const double angleFromHorizontal = d_initAngle_bodyAxis_x + incremented_angle_from_reference_axis[2];
    //double c = cos(angleFromHorizontal), s = sin(angleFromHorizontal);

    // Set the deformation velocity in the body frame.
    //std::vector<double> vec_vel(NDIM);
    //for (std::map<double, int>::const_iterator itr = d_ImmersedBodyData.begin(); itr != d_ImmersedBodyData.end(); itr++)
    //{
    //    double x = itr->first;
    //    const int NumPtsInSection = itr->second;

   for (int lag_idx = 0; lag_idx < d_total_lag_pts; ++lag_idx){
        if (lag_idx < d_large_lag_pts){ // for large sphere
            d_kinematics_vel[0][lag_idx] =  -vc * d_radius_1 * d_radius_1*d_radius_1 / rc3; // x
            d_kinematics_vel[1][lag_idx] =  0.0; //y
            d_kinematics_vel[2][lag_idx] = 0.0;// z
        }else{ // for small sphere
            d_kinematics_vel[0][lag_idx] = vc * d_radius_0 * d_radius_0*d_radius_0 / rc3;;
            d_kinematics_vel[1][lag_idx] = 0.0;;
            d_kinematics_vel[2][lag_idx] = 0.0;//
        }
    }

    return;
} // setEelSpecificVelocity

void
Kinematics::setKinematicsVelocity(const double time,
                                       const std::vector<double>& incremented_angle_from_reference_axis,
                                       const std::vector<double>& center_of_mass,
                                       const std::vector<double>& tagged_pt_position)
{
    d_new_time = time;
    d_incremented_angle_from_reference_axis = incremented_angle_from_reference_axis;
    d_center_of_mass = center_of_mass;
    d_tagged_pt_position = tagged_pt_position;

    setSpecificVelocity(d_new_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);

    return;

} // setNewKinematicsVelocity

const std::vector<std::vector<double> >&
Kinematics::getKinematicsVelocity(const int /*level*/) const
{
    return d_kinematics_vel;

} // getKinematicsVelocity

void
Kinematics::setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis)
{
    return;

    /*printf("Kinematics::setShape time = %g, angle = %g \n", time, incremented_angle_from_reference_axis[2]);

    const StructureParameters& struct_param = getStructureParameters();
    const std::string position_update_method = struct_param.getPositionUpdateMethod();
    */
    // Find the deformed shape. Rotate the shape about center of mass.

    /*TBOX_ASSERT(d_new_time == time);
    double t = time;
    double rc2 = d_radius_0*d_radius_0 + d_radius_1*d_radius_1;
    double new_len = d_amplitude * std::sin(2* PII* d_frequency * t);

    std::vector<double> shape_new(NDIM);

    int lag_idx = -1;

    for (std::map<double, int>::const_iterator itr = d_ImmersedBodyData.begin(); itr != d_ImmersedBodyData.end(); itr++)
    {
        const int NumPtsInSection = itr->second;
        double x = itr->first;
        const double y_shape_base = 0;

        if (x > d_radius_0) {
            x += new_len * d_radius_0 * d_radius_0 / rc2;
        }else {
            x -= new_len * d_radius_1 * d_radius_1 / rc2;
        }

        for (int j = -NumPtsInSection; j <= NumPtsInSection; ++j)
        {
            d_shape[0][++lag_idx] = x;
            d_shape[1][lag_idx] = y_shape_base + j * d_mesh_width[1];
        }

    }
    */
    // Find the c.m of this new shape.
    /*
    std::vector<double> center_of_mass(NDIM, 0.0);
    const int total_lag_pts = d_shape[0].size();
    for (int d = 0; d < NDIM; ++d)
    {
        for (std::vector<double>::const_iterator citr = d_shape[d].begin(); citr != d_shape[d].end(); ++citr)
        {
            center_of_mass[d] += *citr;
        }
    }

    for (int d = 0; d < NDIM; ++d) center_of_mass[d] /= total_lag_pts;

    // Shift the c.m to the origin to apply the rotation
    for (int d = 0; d < NDIM; ++d)
    {
        for (std::vector<double>::iterator itr = d_shape[d].begin(); itr != d_shape[d].end(); ++itr)
        {
            *itr -= center_of_mass[d];
        }
    }
    */
    // Now rotate the shape about origin or center of mass.
    /*
    const double angleFromHorizontal = d_initAngle_bodyAxis_x + d_incremented_angle_from_reference_axis[2];
    for (int i = 0; i < total_lag_pts; ++i)
    {
        const double x_rotated = d_shape[0][i] * cos(angleFromHorizontal) - d_shape[1][i] * sin(angleFromHorizontal);
        const double y_rotated = d_shape[0][i] * sin(angleFromHorizontal) + d_shape[1][i] * cos(angleFromHorizontal);
        d_shape[0][i] = x_rotated;
        d_shape[1][i] = y_rotated;
    }

    d_current_time = d_new_time;

    return;
    */
} // setShape

const std::vector<std::vector<double> >&
Kinematics::getShape(const int /*level*/) const
{
    return d_shape;
} // getShape

} // namespace IBAMR
