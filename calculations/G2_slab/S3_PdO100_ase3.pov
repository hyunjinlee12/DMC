#version 3.6;
#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 2.2 max_trace_level 6}
background {color White}
camera {perspective
  angle 20
  right -9.54*x up 14.60*y
  location <0,0,50.00> look_at <0,0,0>}


light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}
// no fog
#declare simple = finish {phong 0.7 ambient 0.4 diffuse 0.55}
#declare pale = finish {ambient 0.9 diffuse 0.30 roughness 0.001 specular 0.2 }
#declare intermediate = finish {ambient 0.4 diffuse 0.6 specular 0.1 roughness 0.04}
#declare vmd = finish {ambient 0.2 diffuse 0.80 phong 0.25 phong_size 10.0 specular 0.2 roughness 0.1}
#declare jmol = finish {ambient 0.4 diffuse 0.6 specular 1 roughness 0.001 metallic}
#declare ase2 = finish {ambient 0.2 brilliance 3 diffuse 0.6 metallic specular 0.7 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient 0.4 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}
#declare glass = finish {ambient 0.4 diffuse 0.35 specular 1.0 roughness 0.001}
#declare glass2 = finish {ambient 0.3 diffuse 0.3 specular 1.0 reflection 0.25 roughness 0.001}
#declare Rcell = 0.050;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
     torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
     translate LOC}
#end

// no cell vertices
atom(<  0.11,   5.73,  -6.70>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #0
atom(< -3.69,  -2.63, -10.72>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #1
atom(<  0.67,  -3.81,  -6.47>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #2
atom(< -0.75,  -2.48, -11.51>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #3
atom(< -1.56,  -4.62,  -3.16>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #4
atom(< -2.98,  -3.30,  -8.20>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #5
atom(<  1.39,  -4.47,  -3.95>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #6
atom(< -0.04,  -3.14,  -8.99>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #7
atom(< -2.22,  -1.00,  -4.91>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #8
atom(< -3.64,   0.32,  -9.96>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #9
atom(<  0.72,  -0.85,  -5.70>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #10
atom(< -0.70,   0.48, -10.75>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #11
atom(< -1.51,  -1.67,  -2.39>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #12
atom(< -2.93,  -0.34,  -7.44>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #13
atom(<  1.44,  -1.51,  -3.18>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #14
atom(<  0.01,  -0.19,  -8.23>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #15
atom(< -2.17,   1.95,  -4.15>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #16
atom(< -3.59,   3.28,  -9.19>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #17
atom(<  0.77,   2.11,  -4.94>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #18
atom(< -0.65,   3.43,  -9.99>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #19
atom(< -1.46,   1.29,  -1.63>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #20
atom(< -2.88,   2.62,  -6.67>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #21
atom(<  1.49,   1.44,  -2.42>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #22
atom(<  0.06,   2.77,  -7.46>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #23
atom(< -2.12,   4.91,  -3.39>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #24
atom(< -3.54,   6.24,  -8.43>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #25
atom(<  0.82,   5.06,  -4.18>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #26
atom(< -0.60,   6.39,  -9.22>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #27
atom(< -1.41,   4.25,  -0.87>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #28
atom(< -2.83,   5.57,  -5.91>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #29
atom(< -2.27,  -3.96,  -5.68>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #30
atom(<  1.54,   4.40,  -1.66>, 0.56, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #31
atom(< -1.94,  -5.77,  -4.80>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #32
atom(<  0.42,   3.99,   0.00>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #33
atom(< -3.36,  -4.44,  -9.84>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #34
atom(<  1.00,  -5.61,  -5.59>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #35
atom(< -0.42,  -4.29, -10.63>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #36
atom(< -1.15,  -3.55,  -7.33>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #37
atom(<  0.27,  -4.88,  -2.29>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #38
atom(<  1.79,  -3.40,  -8.12>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #39
atom(<  3.21,  -4.72,  -3.08>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #40
atom(< -1.89,  -2.81,  -4.04>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #41
atom(< -3.31,  -1.49,  -9.08>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #42
atom(<  1.05,  -2.66,  -4.83>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #43
atom(< -0.37,  -1.33,  -9.87>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #44
atom(<  0.32,  -1.92,  -1.53>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #45
atom(< -1.10,  -0.59,  -6.57>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #46
atom(<  1.94,   5.47,  -5.83>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #47
atom(<  3.26,  -1.77,  -2.32>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #48
atom(< -1.84,   0.14,  -3.27>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #49
atom(< -3.26,   1.47,  -8.32>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #50
atom(<  1.10,   0.30,  -4.06>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #51
atom(< -0.32,   1.62,  -9.11>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #52
atom(<  0.37,   1.04,  -0.76>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #53
atom(< -1.05,   2.36,  -5.81>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #54
atom(<  3.31,   1.19,  -1.55>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #55
atom(<  1.89,   2.52,  -6.60>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #56
atom(< -1.79,   3.10,  -2.51>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #57
atom(< -3.21,   4.43,  -7.55>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #58
atom(<  1.15,   3.25,  -3.30>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #59
atom(< -0.27,   4.58,  -8.34>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #60
atom(< -1.00,   5.32,  -5.04>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #61
atom(<  1.84,  -0.44,  -7.36>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #62
atom(<  3.36,   4.15,  -0.79>, 1.18, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #63

// no constraints
