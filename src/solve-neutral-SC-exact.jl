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

function solveNeutralSCExact(output_buf, Q, J, theta, v2roc2, v2thetaoc2, v2phioc2, scan_var, scan_mode, target_time_ratio)
#
# This function uses the time-neutral external Schwarzschild polynomials from https://vixra.org/abs/2402.0117
#

  # These are zero for Schwarzschild metric
  global Qse_ratio=bf0
  global J_ratio=bf0

  # spherical velocity mixing
  global sintheta=sin(theta)
  global sin2theta=sintheta * sintheta
  global v2omegaoc2=v2thetaoc2 + (sin2theta * v2phioc2)

  if scan_mode == 0
#
# M(r)
#
    global r=scan_var
    global r2=r * r
    global r3=r2 * r
    global A=bf8 * ref_G3 * r / ref_c6
    global B=-(bf1 - v2omegaoc2) * bf4 * ref_G2 * r2 / ref_c4
    global B2=B * B
    global C=-(v2roc2 + v2omegaoc2) * bf2 * ref_G * r3 / ref_c2
    global discriminant=B2 - (bf4 * A * C)

#@printf("A: %+.3e, B: %+.3e, C: %+.3e, B2: %+.3e, discriminant: %+.3e\n", A, B, C, B2, discriminant)

    # check if A (and then B) is zero
    if abs(A) < BigFloat("1.0E-300") # threshold must be appropriate for the value of setPrecision()
      if abs(B) < BigFloat("1.0E-300") # threshold must be appropriate for the value of setPrecision()
        @printf("Warning, A=B=0, skipping\n")
      else
        @printf("Warning, A=0, using quadratic part of cubic equation\n")
        global A=B
        global B=C
        global B2=B * B
        global C=bf0
        global discriminant=B2 - (bf4 * A * C)
#@printf("A: %+.3e, B: %+.3e, C: %+.3e, B2: %+.3e, discriminant: %+.3e\n", A, B, C, B2, discriminant)
      end # if B
    end #if A

    # check discriminant for type of roots
    if discriminant < BigFloat("-1.0E-300") # threshold must be appropriate for the value of setPrecision()

      # two distinct complex roots
      @printf("Warning: negative discriminant: %+.3e\n, skipping", discriminant)
    elseif abs(discriminant) < BigFloat("1.0E-300") # threshold must be appropriate for the value of setPrecision()

      # one double real root
      global M=(-B + sqrt(discriminant)) / (bf2 * A)
#@printf("double real root: %+.3e\n", r)
      if (M >= bf0) && (M >= M_start) && (M < M_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, pos_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      elseif (M < bf0) && (-M >= M_start) && (-M < M_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, neg_mass\n", r, -M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end # if M
    else # if discriminant
      # two distinct real roots

      # root 1
      global M=(-B + sqrt(discriminant)) / (bf2 * A)
#@printf("root1: %+.3e\n", M)
      if (M >= bf0) && (M >= M_start) && (M < M_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, pos_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      elseif (M < bf0) && (-M >= M_start) && (-M < M_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, neg_mass\n", r, -M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end # if r

      # root2
      global M=(-B - sqrt(discriminant)) / (bf2 * A)
#@printf("root2: %+.3e\n", M)
      if (M >= bf0) && (M >= r_start) && (M < M_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, pos_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      elseif (M < bf0) && (-M >= M_start) && (-M < M_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, neg_mass\n", r, -M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end # if r
    end # if discriminant
  else # if scan_mode

#
# r(M)
#
    global M=scan_var
    global rs=bf2 * ref_G * M / ref_c2
    global rs2=rs * rs
    global rs3=rs2 * rs
    global A=(v2roc2 + v2omegaoc2) * rs
    global B=(bf1 - v2omegaoc2) * rs2
    global B2=B * B
    global C=-rs3
    global discriminant=B2 - (bf4 * A * C)

#@printf("A: %+.3e, B: %+.3e, C: %+.3e, B2: %+.3e, discriminant: %+.3e\n", A, B, C, B2, discriminant)


    # check if A (and then B) is zero
    if abs(A < BigFloat("1.0E-300")) # threshold must be appropriate for the value of setPrecision()
      if abs(B < BigFloat("1.0E-300")) # threshold must be appropriate for the value of setPrecision()
        @printf("Warning, A=B=0, skipping\n")
      else
        @printf("Warning, A=0, using quadratic part of cubic equation\n")
        global A=B
        global B=C
        global B2=B * B
        global C=bf0
        global discriminant=B2 - (bf4 * A * C)
#@printf("A: %+.3e, B: %+.3e, C: %+.3e, B2: %+.3e, discriminant: %+.3e\n", A, B, C, B2, discriminant)
      end # if B
    end #if A
  
    # check discriminant for type of roots
    if discriminant < BigFloat("-1.0E-300") # threshold must be appropriate for the value of setPrecision()

      # two distinct complex roots
      @printf("Warning: negative discriminant: %+.3e, skipping\n", discriminant)
    elseif abs(discriminant) < BigFloat("1.0E-300") # threshold must be appropriate for the value of setPrecision()

      # one double real root
      global r=-B / (bf2 * A)
#@printf("double real root: %+.3e\n", r)
      if (r >= bf0) && (r >= r_start) && (r < r_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, pos_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      elseif (r < bf0) && (-r >= r_start) && (-r < r_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, neg_mass\n", -r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end # if r
    else # if discriminant
      # two distinct real roots

      # root 1
      global r=(-B + sqrt(discriminant)) / (bf2 * A)
#@printf("root1: %+.3e\n", r)
      if (r >= bf0) && (r >= r_start) && (r < r_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, pos_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      elseif (r < bf0) && (-r >= r_start) && (-r < r_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, neg_mass\n", -r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end # if r

      # root2
      global r=(-B - sqrt(discriminant)) / (bf2 * A)
#@printf("root2: %+.3e\n", r)
      if (r >= bf0) && (r >= r_start) && (r < r_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, pos_mass\n", r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      elseif (r < bf0) && (-r >= r_start) && (-r < r_end)
        @printf(output_buf, "%+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, %+.15e, neg_mass\n", -r, M, Q, J, v2roc2, v2thetaoc2, v2phioc2, Qse_ratio, J_ratio)
      end # if r
    end # if discriminant
  end # if scan_mode

end # function solveNeutralSCExact
