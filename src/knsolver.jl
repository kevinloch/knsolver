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

version=1.0
using Printf
setprecision(128)

#
# frequently used integers encoded as BigFloat
#
bf0=BigFloat("0.0")
bf01=BigFloat("0.1")
bf025=BigFloat("0.25")
bf05=BigFloat("0.5")
bf075=BigFloat("0.75")
bf09=BigFloat("0.9")
bf1=BigFloat("1.0")
bf2=BigFloat("2.0")
bf3=BigFloat("3.0")
bf4=BigFloat("4.0")
bf5=BigFloat("5.0")
bf6=BigFloat("6.0")
bf8=BigFloat("8.0")
bf9=BigFloat("9.0")
bf12=BigFloat("12.0")
bf27=BigFloat("27.0")
ref_pi=BigFloat("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870")

#
# reference values from CODATA 2018 unless otherwise noted
#
ref_c=BigFloat("299792458.0")
ref_c2=ref_c * ref_c
inv_ref_c2=bf1 / ref_c2
ref_c3=ref_c2 * ref_c
ref_c4=ref_c3 * ref_c
ref_c5=ref_c4 * ref_c
ref_c6=ref_c4 * ref_c2
ref_G=BigFloat("6.67430E-11")
ref_G2=ref_G * ref_G
ref_G3=ref_G2 * ref_G
ref_e=BigFloat("1.602176634E-19")
ref_e0=BigFloat("8.8541878128E-12")
ref_h=BigFloat("6.62607015E-34")
ref_hbar=ref_h / (bf2 * ref_pi)
ref_hbaro2=ref_hbar / bf2
ref_lP=sqrt(ref_hbar * ref_G / ref_c3)
ref_mP=sqrt(ref_hbar * ref_c / ref_G)
ref_tP=sqrt(ref_hbar * ref_G / ref_c5)
ref_qP=sqrt(bf4 * ref_pi * ref_e0 * ref_hbar * ref_c)
ref_a=ref_e * ref_e / (ref_qP * ref_qP)
ref_golden=(bf1 + sqrt(bf5)) / bf2

#
# shortcuts derived exclusively from reference values
#
dt=ref_tP
c_dt=ref_c * dt
dt2=dt * dt
c2_dt2=ref_c2 * dt2
inv_c2_dt2=bf1 / c2_dt2

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
    global neg_mass_metric_last=target_time_ratio
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

      # solve neutral Kerr-Newman metric (dt = Planck time, all terms divided by c2dt2)
      # positive mass
      pos_mass_metric=(tterms * deltapm) - (rterms / deltapm) - thetaterms - phiterms
      # negative mass (produces exactly identical roots as positive mass but with inverse sign)
      neg_mass_metric=(tterms * deltanm) - (rterms / deltanm) - thetaterms - phiterms

      # detect/display neutral crossing events
      if ((pos_mass_metric > target_time_ratio) && (pos_mass_metric_last < target_time_ratio)) || ((pos_mass_metric < target_time_ratio) && (pos_mass_metric_last > target_time_ratio))
        Qse_ratio=Qse / M
        if J == bf0
          J_ratio=bf0
        else
          J_ratio=abs(M * r * vphi / ref_hbaro2)
        end
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, pos_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end
      if ((neg_mass_metric > target_time_ratio) && (neg_mass_metric_last < target_time_ratio)) || ((neg_mass_metric < target_time_ratio) && (neg_mass_metric_last > target_time_ratio))
        Qse_ratio=Qse / M
        if J == bf0
          J_ratio=bf0
        else
          J_ratio=abs(M * r * vphi / ref_hbaro2)
        end
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, neg_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end

      # update loop vars
      global pos_mass_metric_last=pos_mass_metric
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
    global neg_mass_metric_last=target_time_ratio
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

      # solve neutral Kerr-Newman metric (dt = Planck time, all terms divided by c2dt2)
      # positive mass
      pos_mass_metric=(tterms * deltapm) - (rterms / deltapm) - thetaterms - phiterms
      # negative mass (produces exactly identical roots as positive mass but with inverse sign)
      neg_mass_metric=(tterms * deltanm) - (rterms / deltanm) - thetaterms - phiterms

      # detect/display neutral crossing events
      if ((pos_mass_metric > target_time_ratio) && (pos_mass_metric_last < target_time_ratio)) || ((pos_mass_metric < target_time_ratio) && (pos_mass_metric_last > target_time_ratio))
        Qse_ratio=Qse / M
        if J == bf0
          J_ratio=bf0
        else
          J_ratio=abs(M * r * vphi / ref_hbaro2)
        end
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, pos_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end
      if ((neg_mass_metric > target_time_ratio) && (neg_mass_metric_last < target_time_ratio)) || ((neg_mass_metric < target_time_ratio) && (neg_mass_metric_last > target_time_ratio))
        Qse_ratio=Qse / M
        if J == bf0
          J_ratio=bf0
        else
          J_ratio=abs(M * r * vphi / ref_hbaro2)
        end
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, neg_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end

      # update loop vars
      global pos_mass_metric_last=pos_mass_metric
      global neg_mass_metric_last=neg_mass_metric
      global r *= r_multiplier;
    end # while
  end # if scan_mode
end # function solveKNCoarse()

function main()

################################################
#
# Configuration Options
#
  # target time ratio
  #global target_time_ratio=-bf1 # negative neutral metric (???)
  global target_time_ratio=bf1 # neutral metric
  #global target_time_ratio=bf0 # spaghetti

  # electric charge
  #global Q=bf0
  #global Q=ref_qP
  global Q=bf2 * ref_qP

  # angular momentum
  #global J=bf0
  global J=ref_hbar / bf2
  #global J=-ref_hbar / bf2 # contra-rotating
  #global J=ref_hbar
  #global J=-ref_hbar # contra-rotating

  # polar angle (0 = +z, pi/2 = equator, pi = -z)
  #global theta=bf0          # +z
  #global theta=ref_pi / bf6 # 30 degrees from +z
  #global theta=ref_pi / bf4 # 45 degrees from +z
  #global theta=ref_pi / bf3 # 60 degrees from +z
  global theta=ref_pi / bf2  # equator

  # dimension-specific velocity ratios (v^2/c^2)
  global v2roc2=bf1
  global v2thetaoc2=bf0
  global v2phioc2=bf1

  # r scan range
  global r_start=BigFloat("1.0E-37")
  global r_end=BigFloat("1.0E-33")

  # abs(M) scan range
  global M_start=BigFloat("1.0E-9")
  global M_end=BigFloat("1.0E-5")
################################################

#
# Initialization
#
  # truncate and open output file
  output_file=open("./tmp.csv", "w+")

  # create output buffer
  global output_buf=IOBuffer()

#
# Solve M(r)
#
  global r_increment=BigFloat("1.0E-2") # this precision determines how many sample points per magnitude on the plot
  global r_multiplier=bf1 + r_increment
  global scan_mode=0 # M(r)
  global scan_var=r_start
  @printf("Begin M(r)\n")
  while scan_var < r_end
    @printf("r: %+.3e\n", scan_var)
    solveKNCoarse(output_buf, Q, J, theta, v2roc2, v2thetaoc2, v2phioc2, scan_var, scan_mode, target_time_ratio)
    global scan_var *= r_multiplier;
  end # while
  @printf("End   M(r)\n")
  # write output_buf to output file
  write(output_file, take!(output_buf)) # note: this also clears output_buf so it can be re-used

#
# Solve r(M)
#
  global M_increment=BigFloat("1.0E-2") # this precision determines how many sample points per magnitude on the plot
  global M_multiplier=bf1 + M_increment
  global scan_mode=1 # r(M)
  global scan_var=M_start
  @printf("Begin r(M)\n")
  while scan_var < M_end
    @printf("M: %+.3e\n", scan_var)
    solveKNCoarse(output_buf, Q, J, theta, v2roc2, v2thetaoc2, v2phioc2, scan_var, scan_mode, target_time_ratio)
    global scan_var *= M_multiplier;
  end # while
  @printf("End   r(M)\n")
  # append output_buf to output file
  write(output_file, take!(output_buf))

  # flush and close output file
  close(output_file)

#
# generate plot, and display
#
  # create plot
  #global command=`gnuplot -c ./plotkn.cfg`
  global command=`gnuplot -c ./plotkn-noratios.cfg`
  run(command)

  # display with Preview app`
  global command=`open -a Preview ./tmp.pdf`
  run(command)

end # function main()
