/**
BSD 3-Clause License

Copyright (c) 2018, Vladyslav Usenko and Nikolaus Demmel.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include <bitset>
#include <set>

#include <Eigen/Dense>
#include <sophus/se3.hpp>

#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/relative_pose/CentralRelativePoseSacProblem.hpp>

#include <visnav/camera_models.h>
#include <visnav/common_types.h>

namespace visnav {

void computeEssential(const Sophus::SE3d& T_0_1, Eigen::Matrix3d& E) {
  const Eigen::Vector3d t_0_1 = T_0_1.translation();
  const Eigen::Matrix3d R_0_1 = T_0_1.rotationMatrix();

  // TODO SHEET 3: compute essential matrix
  // UNUSED(E);
  // UNUSED(t_0_1);
  // UNUSED(R_0_1);
  Eigen::Vector3d t_0_1_norm = t_0_1.normalized();
  Eigen::Matrix<double, 3, 3> T_hat;
  T_hat << 0.0f, -t_0_1_norm(2, 0), t_0_1_norm(1, 0), t_0_1_norm(2, 0), 0.0f,
      -t_0_1_norm(0, 0), -t_0_1_norm(1, 0), t_0_1_norm(0, 0), 0.0f;
  E = T_hat * R_0_1;
}

void findInliersEssential(const KeypointsData& kd1, const KeypointsData& kd2,
                          const std::shared_ptr<AbstractCamera<double>>& cam1,
                          const std::shared_ptr<AbstractCamera<double>>& cam2,
                          const Eigen::Matrix3d& E,
                          double epipolar_error_threshold, MatchData& md) {
  md.inliers.clear();

  for (size_t j = 0; j < md.matches.size(); j++) {
    const Eigen::Vector2d p0_2d = kd1.corners[md.matches[j].first];
    const Eigen::Vector2d p1_2d = kd2.corners[md.matches[j].second];

    // TODO SHEET 3: determine inliers and store in md.inliers
    UNUSED(cam1);
    UNUSED(cam2);
    UNUSED(E);
    UNUSED(epipolar_error_threshold);
    UNUSED(p0_2d);
    UNUSED(p1_2d);

    if (abs(cam1->unproject(p0_2d).transpose() * E * cam2->unproject(p1_2d)) <
        epipolar_error_threshold) {
      std::pair<FeatureId, FeatureId> inlier(md.matches[j].first,
                                             md.matches[j].second);
      md.inliers.push_back(inlier);
    }

    std::vector<std::pair<FeatureId, FeatureId>> inliers;
  }
}

void findInliersRansac(const KeypointsData& kd1, const KeypointsData& kd2,
                       const std::shared_ptr<AbstractCamera<double>>& cam1,
                       const std::shared_ptr<AbstractCamera<double>>& cam2,
                       const double ransac_thresh, const int ransac_min_inliers,
                       MatchData& md) {
  md.inliers.clear();
  md.T_i_j = Sophus::SE3d();

  // TODO SHEET 3: Run RANSAC with using opengv's CentralRelativePose and store
  // the final inlier indices in md.inliers and the final relative pose in
  // md.T_i_j (normalize translation). If the number of inliers is smaller than
  // ransac_min_inliers, leave md.inliers empty. Note that if the initial RANSAC
  // was successful, you should do non-linear refinement of the model parameters
  // using all inliers, and then re-estimate the inlier set with the refined
  // model parameters.
  opengv::bearingVectors_t p0_3d_vec;
  opengv::bearingVectors_t p1_3d_vec;
  for (size_t j = 0; j < md.matches.size(); j++) {
    p0_3d_vec.push_back(cam1->unproject(kd1.corners[md.matches[j].first]));
    p1_3d_vec.push_back(cam2->unproject(kd2.corners[md.matches[j].second]));
  }

  // create the central relative adapter
  opengv::relative_pose::CentralRelativeAdapter adapter(p0_3d_vec, p1_3d_vec);
  opengv::transformation_t best_transformation;
  // create a RANSAC object
  opengv::sac::Ransac<
      opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem>
      ransac;
  // create a CentralRelativePoseSacProblem
  std::shared_ptr<
      opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem>
      relposeproblem_ptr(
          new opengv::sac_problems::relative_pose::
              CentralRelativePoseSacProblem(
                  adapter, opengv::sac_problems::relative_pose::
                               CentralRelativePoseSacProblem::NISTER));
  // run ransac
  int maxIterations = 50;
  ransac.sac_model_ = relposeproblem_ptr;
  ransac.threshold_ = ransac_thresh;
  ransac.max_iterations_ = maxIterations;
  ransac.computeModel();
  // get the result
  best_transformation = ransac.model_coefficients_;
  int inlier_num = ransac.inliers_.size();

  // cout << "inlier_num : " << inlier_num << endl;

  // Succeed and then Refine
  // non-linear optimization (using all available correspondences)
  /*
  opengv::rotation_t initial_rotation = best_transformation.block(0, 0, 3, 3);
  opengv::translation_t initial_translation = best_transformation.block(0, 3, 3,
  1); adapter.sett12(initial_translation); adapter.setR12(initial_rotation);
  */
  best_transformation =
      opengv::relative_pose::optimize_nonlinear(adapter, ransac.inliers_);

  inlier_num = ransac.inliers_.size();
  // cout << "optimize inlier_num" << inlier_num << endl;

  // TO Return
  // typedef Eigen::Matrix<double,3,4> transformation_t;
  Eigen::Matrix3d final_rotation = best_transformation.block(0, 0, 3, 3);
  Eigen::Matrix<double, 3, 1> final_translation =
      best_transformation.block(0, 3, 3, 1);
  // cout << "final_translation: " << final_translation << endl;
  final_translation = final_translation.normalized();
  // cout << "final_translation: " << final_translation << endl;
  Sophus::SE3d matrix(final_rotation, final_translation);
  md.T_i_j = matrix;
  // cout << "md.T_i_j : " << md.T_i_j.matrix() << endl;
  // md.T_i_j.matrix3x4() = best_transformation;
  std::vector<int> inliers;
  ransac.sac_model_->selectWithinDistance(ransac.model_coefficients_,
                                          ransac_thresh, inliers);

  // cout << "inliers.size: " << inliers.size();
  if (inliers.size() > ransac_min_inliers) {
    for (size_t i = 0; i < inliers.size(); i++) {
      std::pair<FeatureId, FeatureId> inlier(md.matches[inliers[i]].first,
                                             md.matches[inliers[i]].second);
      md.inliers.push_back(inlier);
    }
  }
}

}  // namespace visnav
