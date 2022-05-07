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

#include <sophus/se3.hpp>

#include <visnav/common_types.h>
using namespace std;
namespace visnav {

double eps = numeric_limits<double>::epsilon();

template <class T>
Eigen::Matrix<T, 3, 3> w_to_w_hat(const Eigen::Matrix<T, 3, 1>& w) {
  Eigen::Matrix<T, 3, 3> w_hat;
  w_hat << 0.0f, -w(2, 0), w(1, 0), w(2, 0), 0.0f, -w(0, 0), -w(1, 0), w(0, 0),
      0.0f;
  return w_hat;
}

template <class T>
// Implement exp for SO(3)
Eigen::Matrix<T, 3, 3> user_implemented_expmap(
    const Eigen::Matrix<T, 3, 1>& xi) {
  // TODO SHEET 1: implement
  // UNUSED(xi);
  Eigen::Matrix<T, 3, 3> rotation;
  Eigen::Matrix<T, 3, 3> I = Eigen::Matrix<T, 3, 3>::Identity();
  Eigen::Matrix<T, 3, 3> w_hat;
  Eigen::Matrix<T, 3, 3> test;
  using std::sqrt;
  T theta = xi.norm();

  w_hat = w_to_w_hat(xi);

  if (theta < eps) {
    rotation = I + (1 - pow(theta, 2) / 6 + pow(theta, 4) / 120) * w_hat +
               (0.5 - pow(theta, 2) / 24 + pow(theta, 4) / 720) * w_hat * w_hat;
  } else {
    rotation = I + w_hat * sin(theta) / theta +
               w_hat * w_hat * (1 - cos(theta)) / (theta * theta);
  }
  return rotation;
}

// Implement log for SO(3)
template <class T>
Eigen::Matrix<T, 3, 1> user_implemented_logmap(
    const Eigen::Matrix<T, 3, 3>& mat) {
  // TODO SHEET 1: implement
  // UNUSED(mat);
  Eigen::Matrix<T, 3, 1> tmp;
  Eigen::Matrix<T, 3, 1> w;
  tmp << mat(2, 1) - mat(1, 2), mat(0, 2) - mat(2, 0), mat(1, 0) - mat(0, 1);
  T trace = mat.trace();
  T theta = acos((trace - 1) * 0.5);
  if (theta < eps) {
    w = 0.5 * tmp * (1 + pow(theta, 2) / 6 + pow(theta, 4) * 7 / 360);
  } else {
    w = 0.5 * tmp * theta / sin(theta);
  }
  return w;
}

// Implement exp for SE(3)
template <class T>
Eigen::Matrix<T, 4, 4> user_implemented_expmap(
    const Eigen::Matrix<T, 6, 1>& xi) {
  // TODO SHEET 1: implement
  Eigen::Matrix<T, 4, 4> rotation;
  Eigen::Matrix<T, 3, 1> v;
  Eigen::Matrix<T, 3, 1> w;
  Eigen::Matrix<T, 3, 3> w_hat;
  Eigen::Matrix<T, 3, 3> I = Eigen::Matrix<T, 3, 3>::Identity();

  // v << xi(0,0), xi(1,0), xi(2,0);
  // w << xi(3,0), xi(4,0), xi(5,0);
  v = xi.block(0, 0, 3, 1);
  w = xi.block(3, 0, 3, 1);
  w_hat = w_to_w_hat(w);

  T theta = w.norm();

  Eigen::Matrix<T, 3, 3> exp_w_hat = user_implemented_expmap(w);
  Eigen::Matrix<T, 3, 3> J;
  if (theta < eps) {
    J = I + (1 - (0.5 - pow(theta, 2) / 24 + pow(theta, 4) / 720)) * w_hat +
        (1 / 6 - pow(theta, 2) / 120 + pow(theta, 4) / 5040) * w_hat * w_hat;
  } else {
    J = I + ((1 - cos(theta)) / (theta * theta)) * w_hat +
        ((theta - sin(theta)) / (theta * theta * theta)) * w_hat * w_hat;
  }
  Eigen::Matrix<T, 3, 1> Jv = J * v;

  /*
  rotation << exp_w_hat(0,0), exp_w_hat(0,1), exp_w_hat(0,2), Jv(0,0),
              exp_w_hat(1,0), exp_w_hat(1,1), exp_w_hat(1,2), Jv(1,0),
              exp_w_hat(2,0), exp_w_hat(2,1), exp_w_hat(2,2), Jv(2,0),
              0,0,0,1;
  */
  rotation.block(0, 0, 3, 3) = exp_w_hat;
  rotation.block(0, 3, 3, 1) = Jv;
  rotation.block(3, 0, 1, 4) << 0, 0, 0, 1;

  return rotation;
}

// Implement log for SE(3)
template <class T>
Eigen::Matrix<T, 6, 1> user_implemented_logmap(
    const Eigen::Matrix<T, 4, 4>& mat) {
  // TODO SHEET 1: implement
  // UNUSED(mat);
  Eigen::Matrix<T, 3, 3> R;
  Eigen::Matrix<T, 3, 3> I = Eigen::Matrix<T, 3, 3>::Identity();
  Eigen::Matrix<T, 3, 1> v;
  Eigen::Matrix<T, 6, 1> twist;
  Eigen::Matrix<T, 3, 3> J_inverse;

  /*
  R << mat(0,0), mat(0,1), mat(0,2),
       mat(1,0), mat(1,1), mat(1,2),
       mat(2,0), mat(2,1), mat(2,2);
  */
  R = mat.block(0, 0, 3, 3);
  Eigen::Matrix<T, 3, 1> w = user_implemented_logmap(R);
  Eigen::Matrix<T, 3, 1> t;
  // t << mat(0,3), mat(1,3), mat(2,3);
  t = mat.block(0, 3, 3, 1);

  T theta = w.norm();
  Eigen::Matrix<T, 3, 3> w_hat;
  w_hat = w_to_w_hat(w);

  if (theta < eps) {
    J_inverse =
        I - 0.5 * w_hat +
        (1 / 12 + pow(theta, 2) / 720 + pow(theta, 4) / 30240) * w_hat * w_hat;
  } else {
    J_inverse =
        I - 0.5 * w_hat +
        (1 / pow(theta, 2) - (1 + cos(theta)) / (2 * theta * sin(theta))) *
            w_hat * w_hat;
  }
  v = J_inverse * t;
  twist << v(0, 0), v(1, 0), v(2, 0), w(0, 0), w(1, 0), w(2, 0);
  return twist;
}

}  // namespace visnav
