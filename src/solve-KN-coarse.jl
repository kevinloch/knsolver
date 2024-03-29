#
# knsolver is a tool for exploring/analyzing the Kerr-Newman metric near Planck scale
# Configuration options are in main()
#

# BSD 3-Clause License
#
# Copyright (c) 2024, Kevin Loch
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function solveKNCoarse(output_buf, Q, J, theta, v2roc2, v2thetaoc2, v2phioc2, scan_var, scan_mode, target_time_ratio)
  global rq2=Q * Q * ref_G / (bf4 * ref_pi * ref_e0 * ref_c4)
  global dr2=v2roc2 * ref_c2 * dt2
  global sintheta=sin(theta)
  global costheta=cos(theta)
  global sin2theta=sintheta * sintheta
  global cos2theta=costheta * costheta
  vphi=sqrt(v2phioc2 * ref_c2)

  if scan_mode == 0
#
# M(r)
#
    global r=scan_var
    global Qse=Q * Q / (bf4 * ref_pi * ref_e0 * ref_c2 * r) # Electromagnetic self-energy mass equivalent
    global r2=r * r
    global inv_r2 = bf1 / r2
    global dtheta2=v2thetaoc2 * ref_c2 * dt2 * inv_r2
    global dphi=vphi * dt / r

    # inner loop to scan M
    global M=M_start
    global M_increment=BigFloat("1.0E-3") # this precision determines how close together different roots can be resolved
    global M_multiplier=bf1 + M_increment
    global pos_mass_metric_last=target_time_ratio
    global pos_mass_metric_nlast=target_time_ratio
    global neg_mass_metric_last=target_time_ratio
    global neg_mass_metric_nlast=target_time_ratio
    while M <= M_end
      rs=bf2 * ref_G * M * inv_ref_c2
      rs_r=rs * r
      a=J / (M * ref_c)
      a2=a * a
      rho2=r2 + (a2 * cos2theta)
      inv_rho2_c2_dt2=bf1 / (rho2 * c2_dt2)
      c_dt_m_a_sin2theta_dphi=c_dt - (a * sin2theta * dphi)
      tterms=c_dt_m_a_sin2theta_dphi * c_dt_m_a_sin2theta_dphi  * inv_rho2_c2_dt2
      rterms=dr2 * rho2 * inv_c2_dt2
      thetaterms=dtheta2 * rho2 * inv_c2_dt2
      r2_p_a2_dphi_m_a_c_dt=((r2 + a2) * dphi) - (a * c_dt)
      phiterms=r2_p_a2_dphi_m_a_c_dt * r2_p_a2_dphi_m_a_c_dt * sin2theta * inv_rho2_c2_dt2
      deltapm=r2 - rs_r + a2 + rq2
      deltanm=r2 + rs_r + a2 + rq2

      # solve Kerr-Newman metric (dt = Planck time, all terms divided by c2dt2)
      # positive mass
      pos_mass_metric=(tterms * deltapm) - (rterms / deltapm) - thetaterms - phiterms
      # negative mass (produces exactly identical roots as positive mass but with inverse sign)
      neg_mass_metric=(tterms * deltanm) - (rterms / deltanm) - thetaterms - phiterms

      # detect/display target time ratio crossing events, but ignore divergent crossings by tracking previous slope with *nlast
      if ((pos_mass_metric > target_time_ratio) && (pos_mass_metric_last < target_time_ratio) && (pos_mass_metric_nlast < pos_mass_metric_last)) || ((pos_mass_metric < target_time_ratio) && (pos_mass_metric_last > target_time_ratio) && (pos_mass_metric_nlast > pos_mass_metric_last))
        Qse_ratio=Qse / M
        if J == bf0
          J_ratio=bf0
        else
          J_ratio=abs(M * r * vphi / ref_hbaro2)
        end
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, pos_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end
      # detect/display target time ratio crossing events, but ignore divergent crossings by tracking previous slope with *nlast
      if ((neg_mass_metric > target_time_ratio) && (neg_mass_metric_last < target_time_ratio) && (neg_mass_metric_nlast < neg_mass_metric_last)) || ((neg_mass_metric < target_time_ratio) && (neg_mass_metric_last > target_time_ratio) && (neg_mass_metric_nlast > neg_mass_metric_last))
        Qse_ratio=Qse / M
        if J == bf0
          J_ratio=bf0
        else
          J_ratio=abs(M * r * vphi / ref_hbaro2)
        end
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, neg_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end

      # update loop vars
      global pos_mass_metric_nlast=pos_mass_metric_last
      global pos_mass_metric_last=pos_mass_metric
      global neg_mass_metric_nlast=neg_mass_metric_last
      global neg_mass_metric_last=neg_mass_metric
      global M *= M_multiplier;
    end # while
  else # if scan_mode
#
# r(M)
#
    global M=scan_var
    global rs=bf2 * ref_G * M * inv_ref_c2
    global a=J / (M * ref_c)
    global a2=a * a
    global a_c_dt=a * c_dt

    # inner loop to scan r
    global r=r_start
    global r_increment=BigFloat("1.0E-3") # this precision determines how close together different roots can be resolved
    global r_multiplier=bf1 + r_increment
    global pos_mass_metric_last=target_time_ratio
    global pos_mass_metric_nlast=target_time_ratio
    global neg_mass_metric_last=target_time_ratio
    global neg_mass_metric_nlast=target_time_ratio
    while r <= r_end
      Qse=Q * Q / (bf4 * ref_pi * ref_e0 * ref_c2 * r) # Electromagnetic self-energy mass equivalent
      r2=r * r
      inv_r2 = bf1 / r2
      dtheta2=v2thetaoc2 * ref_c2 * dt2 * inv_r2
      dphi=vphi * dt / r
      rs_r=rs * r
      rho2=r2 + (a2 * cos2theta)
      inv_rho2_c2_dt2=bf1 / (rho2 * c2_dt2)
      c_dt_m_a_sin2theta_dphi=c_dt - (a * sin2theta * dphi)
      tterms=c_dt_m_a_sin2theta_dphi * c_dt_m_a_sin2theta_dphi  * inv_rho2_c2_dt2
      rterms=dr2 * rho2 * inv_c2_dt2
      thetaterms=dtheta2 * rho2 * inv_c2_dt2
      r2_p_a2_dphi_m_a_c_dt=((r2 + a2) * dphi) - a_c_dt
      phiterms=r2_p_a2_dphi_m_a_c_dt * r2_p_a2_dphi_m_a_c_dt * sin2theta * inv_rho2_c2_dt2
      deltapm=r2 - rs_r + a2 + rq2
      deltanm=r2 + rs_r + a2 + rq2

      # solve Kerr-Newman metric (dt = Planck time, all terms divided by c2dt2)
      # positive mass
      pos_mass_metric=(tterms * deltapm) - (rterms / deltapm) - thetaterms - phiterms
      # negative mass (produces exactly identical roots as positive mass but with inverse sign)
      neg_mass_metric=(tterms * deltanm) - (rterms / deltanm) - thetaterms - phiterms

      # detect/display target time ratio crossing events, but ignore divergent crossings by tracking previous slope with *nlast
      if ((pos_mass_metric > target_time_ratio) && (pos_mass_metric_last < target_time_ratio) && (pos_mass_metric_nlast < pos_mass_metric_last)) || ((pos_mass_metric < target_time_ratio) && (pos_mass_metric_last > target_time_ratio) && (pos_mass_metric_nlast > pos_mass_metric_last))
        Qse_ratio=Qse / M
        if J == bf0
          J_ratio=bf0
        else
          J_ratio=abs(M * r * vphi / ref_hbaro2)
        end
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, pos_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end
      # detect/display target time ratio crossing events, but ignore divergent crossings by tracking previous slope with *nlast
      if ((neg_mass_metric > target_time_ratio) && (neg_mass_metric_last < target_time_ratio) && (neg_mass_metric_nlast < neg_mass_metric_last)) || ((neg_mass_metric < target_time_ratio) && (neg_mass_metric_last > target_time_ratio) && (neg_mass_metric_nlast > neg_mass_metric_last))
        Qse_ratio=Qse / M
        if J == bf0
          J_ratio=bf0
        else
          J_ratio=abs(M * r * vphi / ref_hbaro2)
        end
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, neg_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end

      # update loop vars
      global pos_mass_metric_nlast=pos_mass_metric_last
      global pos_mass_metric_last=pos_mass_metric
      global neg_mass_metric_nlast=neg_mass_metric_last
      global neg_mass_metric_last=neg_mass_metric
      global r *= r_multiplier;
    end # while
  end # if scan_mode
end # function solveKNCoarse()
