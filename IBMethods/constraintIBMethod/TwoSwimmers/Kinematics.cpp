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

double distance(std::vector<double> &u, std::vector<double> &v, int DIM)
{
    double dist = 0.0;
    for (int d = 0; d < DIM; ++d){
        dist += (u[d] - v[d])*(u[d] - v[d]);
    }
    return std::sqrt(dist);
}

double angle(std::vector<double> &u, std::vector<double> &v, int DIM)
{
    double costheta, uv, umag, vmag;
    uv = umag = vmag = 0.0;

    for (int d = 0; d < DIM; ++d){
        uv   += u[d] * v[d];
        umag += u[d] * u[d];
        vmag += v[d] * v[d];
    }

    costheta = uv / (std::sqrt(umag) * std::sqrt(vmag));
    return std::acos(costheta);
}
} // namespace

///////////////////////////////////////////////////////////////////////

Kinematics::Kinematics(const std::string& object_name,
                                 Pointer<Database> input_db,
                                 LDataManager* l_data_manager,
                                 Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                 double ds,
                                 bool register_for_restart)
    : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
      d_current_time(0.0),
      d_kinematics_vel(NDIM),
      d_shape(NDIM),
      d_center_of_mass(3),
      d_incremented_angle_from_reference_axis(3),
      d_tagged_pt_position(3),
      d_mesh_width(NDIM),
      d_ds(ds)
{
    std::vector<double> center_0(NDIM), center_1(NDIM);
    input_db->getDoubleArray("center_0", &center_0[0], NDIM);
    input_db->getDoubleArray("center_1", &center_1[0], NDIM);

    d_radius_0  = input_db->getDouble("radius_0");
    d_radius_1  = input_db->getDouble("radius_1");
    d_amplitude = input_db->getDouble("amplitude");
    d_frequency = input_db->getDouble("frequency");

    std::vector<double> cc_vec(NDIM), xAxis{1.0, 0, 0};
    for (int d = 0; d < NDIM; ++d){
      cc_vec[d] = center_1[d] - center_0[d];
    }

    d_center_to_center = distance(center_0, center_1, NDIM);
    double iniAngle = angle(cc_vec, xAxis, NDIM);
    if (cc_vec[1] < 0.0) iniAngle = 2*PII - iniAngle;
    d_initAngle_bodyAxis_x = iniAngle;

    printf("Kinematics::Kinematics() ini angle wrt x-axis = %g\n", d_initAngle_bodyAxis_x * 180 /PII);

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
        d_mesh_width[dim] = d_ds * dx[dim];
    }

    const int Disc_0_Nx = static_cast<int>(round(d_radius_0 / d_mesh_width[0]));
    const int Disc_1_Nx = static_cast<int>(round(d_radius_1 / d_mesh_width[0]));

    d_accepted_radius_0 = Disc_0_Nx * d_mesh_width[0];
    d_accepted_radius_1 = Disc_1_Nx * d_mesh_width[0];

    d_ImmersedBodyData.clear();
    d_n_pts_0 = 1;
    d_n_pts_1 = 1;

    // center 0, we encode based on ring radius, for disc 1, we add radius 0 to the pair
    d_ImmersedBodyData.insert(std::make_pair(0, 1));

    for (int i = 1; i <= Disc_0_Nx; ++i)
    {
        double r = i * d_mesh_width[0];
        int nPtsQter = static_cast<int>(round(0.5 * PII * r / d_mesh_width[0]));

        d_ImmersedBodyData.insert(std::make_pair(r, 4*nPtsQter));
        d_n_pts_0 += 4*nPtsQter;
    }

    // center 1 starting from radius_0 + ds
    double r_base = d_accepted_radius_0 + d_mesh_width[0];
    d_ImmersedBodyData.insert(std::make_pair(r_base, 1));

    for (int i = 1; i <= Disc_1_Nx; ++i)
    {
        double r =  i * d_mesh_width[0];
        int nPtsQter = static_cast<int>(round(0.5 * PII * r / d_mesh_width[0]));
        d_ImmersedBodyData.insert(std::make_pair(r + r_base, 4*nPtsQter));
        d_n_pts_1 += 4*nPtsQter;
    }

    printf("Kinematics::setImmersedBodyLayout() lag index range = [%d %d], total lag pts = %d\n",\
     idx_range[0].first, idx_range[0].second, total_lag_pts);
    printf("Kinematics::setImmersedBodyLayout() mesh width = [%g, %g]\n", d_mesh_width[0], d_mesh_width[1]);
    printf("Kinematics::setImmersedBodyLayout() vertex number = [%d, %d]\n", d_n_pts_0, d_n_pts_1);
    printf("Kinematics::setImmersedBodyLayout() accepted radius = [%g, %g]\n", d_accepted_radius_0, d_accepted_radius_1);

    return;

} // setImmersedBodyLayout

void
Kinematics::setSpecificVelocity(const double time,
                                        const std::vector<double>& incremented_angle_from_reference_axis,
                                        const std::vector<double>& center_of_mass,
                                        const std::vector<double>& tagged_pt_position)
{
    //double t  = time;
    double vc = d_amplitude * 2 * PII * d_frequency * std::cos(2* PII* d_frequency * time);
    double rr = d_accepted_radius_0 * d_accepted_radius_0 + d_accepted_radius_1 * d_accepted_radius_1;
    double len_scale_0 = d_accepted_radius_1 * d_accepted_radius_1 / rr;
    double len_scale_1 = d_accepted_radius_0 * d_accepted_radius_0 / rr;

    const double angleFromHorizontal = d_initAngle_bodyAxis_x + incremented_angle_from_reference_axis[2];
    double c = cos(angleFromHorizontal), s = sin(angleFromHorizontal);



    // Set the deformation velocity in the body frame.
    std::vector<double> vec_vel(NDIM);

    int lag_idx = 0, lag_small = 0, lag_large = 0;
    for (std::map<double, int>::const_iterator itr = d_ImmersedBodyData.begin(); itr != d_ImmersedBodyData.end(); itr++)
    {
        const double r = itr->first;
        const int nPts = itr->second;
        // small disc
        if (r > d_accepted_radius_0){
            vec_vel[0] =  c * vc * len_scale_1;
            vec_vel[1] =  s * vc * len_scale_1;
            lag_small += nPts;
        }else{
            vec_vel[0] = -c * vc * len_scale_0;
            vec_vel[1] = -s * vc * len_scale_0;
            lag_large += nPts;
        }

        const int lowerlimit = lag_idx;
        const int upperlimit = lag_idx + nPts;
        for (int d = 0; d < NDIM; ++d)
        {
            for (int i = lowerlimit; i < upperlimit; ++i) d_kinematics_vel[d][i] = vec_vel[d];
        }

        lag_idx = upperlimit;
    }

    //printf("Kinematics::setSpecificVelocity() time %g, current angle wrt x-axis %g (deg)\n", time, angleFromHorizontal * 180/PII);
    //printf("Kinematics::setSpecificVelocity() time %g, lag pts = [%d, %d]\n", time, lag_large, lag_small);
    return;
} // setSpecificVelocity

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

const std::vector<double>&
Kinematics::getCurrentDiscPosition()
{
    Eigen::Vector2d com, r0, r1;
    double angle_from_x = d_initAngle_bodyAxis_x + d_incremented_angle_from_reference_axis[2];
    double bot_len = d_center_to_center + d_amplitude * std::sin(2* PII* d_frequency * d_new_time);
    double len_0, len_1;

    //amr lengh of each disc ~ area of each disc
    //double rr = d_accepted_radius_0*d_accepted_radius_0 + d_accepted_radius_1*d_accepted_radius_1;
    //len_0 =  bot_len * d_accepted_radius_1 * d_accepted_radius_1 / rr;
    //len_1 =  bot_len * d_accepted_radius_0 * d_accepted_radius_0 / rr;

    // amr lengh of each disc ~ number of lag points
    len_0 =  bot_len * d_n_pts_1 /(d_n_pts_0 + d_n_pts_1);
    len_1 =  bot_len * d_n_pts_0 /(d_n_pts_0 + d_n_pts_1);

    r1 <<  len_1 * std::cos(angle_from_x),  len_1 * std::sin(angle_from_x);
    r0 << -len_0 * std::cos(angle_from_x),   -len_0 * std::sin(angle_from_x);

    com << d_center_of_mass[0], d_center_of_mass[1];

    r1 += com;
    r0 += com;

    d_disc_pos.clear();

    d_disc_pos.push_back(r0(0));
    d_disc_pos.push_back(r0(1));
    d_disc_pos.push_back(r1(0));
    d_disc_pos.push_back(r1(1));

    d_disc_pos.push_back(d_n_pts_0);
    d_disc_pos.push_back(d_n_pts_1);

    //printf("points %d, %d, tot = %d\n", d_n_pts_0, d_n_pts_1, d_n_pts_0 + d_n_pts_1);
    /*
    printf("Kinematics::getCurrentDiscPosition() time %g struct com = [%g, %g], bot_len %g\n", d_new_time, com(0), com(1), bot_len);
    printf("Kinematics::getCurrentDiscPosition() time %g disc position: large = [%g, %g], small = [%g, %g]\n", \
    d_new_time, r0(0), r0(1), r1(0), r1(1));*/
    return d_disc_pos;
}

const std::vector<std::vector<double> >&
Kinematics::getKinematicsVelocity(const int /*level*/) const
{
    return d_kinematics_vel;

} // getKinematicsVelocity

void
Kinematics::setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis)
{

    const StructureParameters& struct_param = getStructureParameters();
    const std::string position_update_method = struct_param.getPositionUpdateMethod();
    if (position_update_method == "CONSTRAINT_VELOCITY") return;

    // Find the deformed shape. Rotate the shape about center of mass.
    TBOX_ASSERT(d_new_time == time);
    double t = time;
    double rr = d_accepted_radius_0*d_accepted_radius_0 + d_accepted_radius_1*d_accepted_radius_1;
    double new_len = d_amplitude * std::sin(2* PII* d_frequency * t);
    double len_scale_0 = d_accepted_radius_1 * d_accepted_radius_1 / rr;
    double len_scale_1 = d_accepted_radius_0 * d_accepted_radius_0 / rr;

    int lag_idx = -1, lag_large = 0, lag_small = 0;
    for (std::map<double, int>::const_iterator itr = d_ImmersedBodyData.begin(); itr != d_ImmersedBodyData.end(); itr++)
    {
        const int nPts = itr->second;
        const double r_shifted = itr->first;
        double cent_x, r, disp_x, phi;

        // small disc
        if (r_shifted > d_accepted_radius_0) {
            r      =  r_shifted - d_accepted_radius_0 - d_mesh_width[0];
            disp_x = -new_len * len_scale_1;
            cent_x = d_center_to_center;
            lag_small += nPts;
        // large disc
        }else {
            r      = r_shifted;
            disp_x = new_len * len_scale_0;
            cent_x = 0.0;
            lag_large += nPts;
        }

        for (int j = 0; j < nPts; ++j)
        {
            phi = j * 2*PII/nPts;
            d_shape[0][++lag_idx] = r * std::cos(phi) + cent_x + disp_x;
            d_shape[1][lag_idx]   = r * std::sin(phi);
        }
    }

    // Find the c.m of this new shape.
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

    // Now rotate the shape about origin or center of mass.
    const double angleFromHorizontal = d_initAngle_bodyAxis_x + d_incremented_angle_from_reference_axis[2];
    printf("Kinematics::setShape() time %g, current angle wrt x-axis %g(deg)\n", t, angleFromHorizontal * 180/PII);
    printf("Kinematics::setShape() time %g cc_dist %g, lag pts [%d, %d]\n", t, d_center_to_center, lag_large, lag_small);

    for (int i = 0; i < total_lag_pts; ++i)
    {
        const double x_rotated = d_shape[0][i] * cos(angleFromHorizontal) - d_shape[1][i] * sin(angleFromHorizontal);
        const double y_rotated = d_shape[0][i] * sin(angleFromHorizontal) + d_shape[1][i] * cos(angleFromHorizontal);
        d_shape[0][i] = x_rotated;
        d_shape[1][i] = y_rotated;
    }

    d_current_time = d_new_time;

    return;
} // setShape

const std::vector<std::vector<double> >&
Kinematics::getShape(const int /*level*/) const
{
    return d_shape;
} // getShape

} // namespace IBAMR
