/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     5.0                                |
|   \\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\/     M anipulation  | For copyright notice see file Copyright         |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//process this file using: m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions -----------------------------
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'printf ($1)')])
m4_define(pi, 3.14159265358979323844)
m4_define(rad, [calc($1*pi/180.0)])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

//Geometry -----------------------------------
// 2 planes levels for rotor
m4_define(zA, 0.0)
m4_define(zB, 0.5)

// 2 planes levels for stator
m4_define(zC, 0.0)
m4_define(zD, 0.5)

// Angle span for inner rotor block
m4_define(angleStartRotor, rad(  0.0))
m4_define(angleStopRotor,  rad(360.0))

// Angle span for outer stator block
m4_define(angleStartStator, rad( 0.0))
m4_define(angleStopStator, rad( 360.0))

// Radial dimensions
m4_define(rMinRotor,      1.0)
m4_define(rRotorImpeller, 2.0)  // Impeller: from rMinRotor to rRotorImpeller
m4_define(rMaxRotor,      3.0)

m4_define(rMinStator,     3.0)
m4_define(rStatorBaffle,  4.0)  // Baffle : from rStatorBaffle to rMaxStator
m4_define(rMaxStator,     5.0)

// Baffle/impeller half width angle
m4_define(impThick, rad(0.60))  // impeller thickness
m4_define(bafThick, rad(0.12))  // baffle thickness

// Mesh parameters
m4_define(BLOCKSIZE_ROTOR_SECTOR_090,   8  27 4)
m4_define(BLOCKSIZE_STATOR_SECTOR_060,  8  18 4)
m4_define(grading, 1.0)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.1;

vertices
(
// Stator with impeller
//Plane A:
//Bottom of rotor block#1 : sector 0-90 deg
(calc(rMinRotor     *cos(angleStartRotor+impThick))                  calc(rMinRotor     *sin(angleStartRotor+impThick))                   zA) vlabel(A0_000)
(calc(rRotorImpeller*cos(angleStartRotor))                           calc(rRotorImpeller*sin(angleStartRotor))                            zA) vlabel(A1_000)
(calc(rMinRotor     *cos(angleStartRotor+rad(90)-impThick))          calc(rMinRotor     *sin(angleStartRotor+rad(90)-impThick))           zA) vlabel(A2_000)
(calc(rRotorImpeller*cos(angleStartRotor+rad(90)))                   calc(rRotorImpeller*sin(angleStartRotor+rad(90)))                    zA) vlabel(A3_000)

//Plane B:
//Top of rotor block#1 : sector 0-90 deg
(calc(rMinRotor     *cos(angleStartRotor+impThick))                  calc(rMinRotor     *sin(angleStartRotor+impThick))                   zB) vlabel(B0_000)
(calc(rRotorImpeller*cos(angleStartRotor))                           calc(rRotorImpeller*sin(angleStartRotor))                            zB) vlabel(B1_000)
(calc(rMinRotor     *cos(angleStartRotor+rad(90)-impThick))          calc(rMinRotor     *sin(angleStartRotor+rad(90)-impThick))           zB) vlabel(B2_000)
(calc(rRotorImpeller*cos(angleStartRotor+rad(90)))                   calc(rRotorImpeller*sin(angleStartRotor+rad(90)))                    zB) vlabel(B3_000)

//Plane A:
//Bottom of rotor block#2 : sector 90-180 deg
(calc(rMinRotor     *cos(angleStartRotor+impThick+rad(90)))          calc(rMinRotor     *sin(angleStartRotor+impThick+rad(90)))           zA) vlabel(A0_090)
(calc(rRotorImpeller*cos(angleStartRotor+rad(90)))                   calc(rRotorImpeller*sin(angleStartRotor+rad(90)))                    zA) vlabel(A1_090)
(calc(rMinRotor     *cos(angleStartRotor+rad(90)-impThick+rad(90)))  calc(rMinRotor     *sin(angleStartRotor+rad(90)-impThick+rad(90)))   zA) vlabel(A2_090)
(calc(rRotorImpeller*cos(angleStartRotor+rad(90)+rad(90)))           calc(rRotorImpeller*sin(angleStartRotor+rad(90)+rad(90)))            zA) vlabel(A3_090)

//Plane B:
//Top of rotor block#2 : sector 90-180 deg
(calc(rMinRotor     *cos(angleStartRotor+impThick+rad(90)))          calc(rMinRotor     *sin(angleStartRotor+impThick+rad(90)))           zB) vlabel(B0_090)
(calc(rRotorImpeller*cos(angleStartRotor+rad(90)))                   calc(rRotorImpeller*sin(angleStartRotor+rad(90)))                    zB) vlabel(B1_090)
(calc(rMinRotor     *cos(angleStartRotor+rad(90)-impThick+rad(90)))  calc(rMinRotor     *sin(angleStartRotor+rad(90)-impThick+rad(90)))   zB) vlabel(B2_090)
(calc(rRotorImpeller*cos(angleStartRotor+rad(90)+rad(90)))           calc(rRotorImpeller*sin(angleStartRotor+rad(90)+rad(90)))            zB) vlabel(B3_090)

//Plane A:
//Bottom of rotor block#3 : sector 180-270 deg
(calc(rMinRotor     *cos(angleStartRotor+impThick+rad(180)))         calc(rMinRotor     *sin(angleStartRotor+impThick+rad(180)))          zA) vlabel(A0_180)
(calc(rRotorImpeller*cos(angleStartRotor+rad(180)))                  calc(rRotorImpeller*sin(angleStartRotor+rad(180)))                   zA) vlabel(A1_180)
(calc(rMinRotor     *cos(angleStartRotor+rad(90)-impThick+rad(180))) calc(rMinRotor     *sin(angleStartRotor+rad(90)-impThick+rad(180)))  zA) vlabel(A2_180)
(calc(rRotorImpeller*cos(angleStartRotor+rad(90)+rad(180)))          calc(rRotorImpeller*sin(angleStartRotor+rad(90)+rad(180)))           zA) vlabel(A3_180)

//Plane B:
//Top of rotor block#3 : sector 180-270 deg
(calc(rMinRotor     *cos(angleStartRotor+impThick+rad(180)))         calc(rMinRotor     *sin(angleStartRotor+impThick+rad(180)))          zB) vlabel(B0_180)
(calc(rRotorImpeller*cos(angleStartRotor+rad(180)))                  calc(rRotorImpeller*sin(angleStartRotor+rad(180)))                   zB) vlabel(B1_180)
(calc(rMinRotor     *cos(angleStartRotor+rad(90)-impThick+rad(180))) calc(rMinRotor     *sin(angleStartRotor+rad(90)-impThick+rad(180)))  zB) vlabel(B2_180)
(calc(rRotorImpeller*cos(angleStartRotor+rad(90)+rad(180)))          calc(rRotorImpeller*sin(angleStartRotor+rad(90)+rad(180)))           zB) vlabel(B3_180)

//Plane A:
//Bottom of rotor block#4 : sector 270-360 deg
(calc(rMinRotor     *cos(angleStartRotor+impThick+rad(270)))         calc(rMinRotor     *sin(angleStartRotor+impThick+rad(270)))          zA) vlabel(A0_270)
(calc(rRotorImpeller*cos(angleStartRotor+rad(270)))                  calc(rRotorImpeller*sin(angleStartRotor+rad(270)))                   zA) vlabel(A1_270)
(calc(rMinRotor     *cos(angleStartRotor+rad(90)-impThick+rad(270))) calc(rMinRotor     *sin(angleStartRotor+rad(90)-impThick+rad(270)))  zA) vlabel(A2_270)
(calc(rRotorImpeller*cos(angleStartRotor+rad(90)+rad(270)))          calc(rRotorImpeller*sin(angleStartRotor+rad(90)+rad(270)))           zA) vlabel(A3_270)

//Plane B:
//Top of rotor block#4 : sector 270-360 deg
(calc(rMinRotor     *cos(angleStartRotor+impThick+rad(270)))         calc(rMinRotor     *sin(angleStartRotor+impThick+rad(270)))          zB) vlabel(B0_270)
(calc(rRotorImpeller*cos(angleStartRotor+rad(270)))                  calc(rRotorImpeller*sin(angleStartRotor+rad(270)))                   zB) vlabel(B1_270)
(calc(rMinRotor     *cos(angleStartRotor+rad(90)-impThick+rad(270))) calc(rMinRotor     *sin(angleStartRotor+rad(90)-impThick+rad(270)))  zB) vlabel(B2_270)
(calc(rRotorImpeller*cos(angleStartRotor+rad(90)+rad(270)))          calc(rRotorImpeller*sin(angleStartRotor+rad(90)+rad(270)))           zB) vlabel(B3_270)

//Plane A:
//Bottom of rotor block#5 : sector without impeller 0-90 deg
(calc(rMaxRotor     *cos(angleStartRotor))                           calc(rMaxRotor     *sin(angleStartRotor))                            zA) vlabel(AA1_000)
(calc(rMaxRotor     *cos(angleStartRotor+rad(90)))                   calc(rMaxRotor     *sin(angleStartRotor+rad(90)))                    zA) vlabel(AA3_000)

//Plane B:
//Top of rotor block#5 : sector without impeller 0-90 deg
(calc(rMaxRotor     *cos(angleStartRotor))                           calc(rMaxRotor     *sin(angleStartRotor))                            zB) vlabel(BB1_000)
(calc(rMaxRotor     *cos(angleStartRotor+rad(90)))                   calc(rMaxRotor     *sin(angleStartRotor+rad(90)))                    zB) vlabel(BB3_000)

//Plane A:
//Bottom of rotor block#6 : sector without impeller 90-180 deg
(calc(rMaxRotor     *cos(angleStartRotor+rad(90)))                   calc(rMaxRotor     *sin(angleStartRotor+rad(90)))                    zA) vlabel(AA1_090)
(calc(rMaxRotor     *cos(angleStartRotor+rad(90)+rad(90)))           calc(rMaxRotor     *sin(angleStartRotor+rad(90)+rad(90)))            zA) vlabel(AA3_090)

//Plane B:
//Top of rotor block#6 : sector without impeller 90-180 deg
(calc(rMaxRotor     *cos(angleStartRotor+rad(90)))                   calc(rMaxRotor     *sin(angleStartRotor+rad(90)))                    zB) vlabel(BB1_090)
(calc(rMaxRotor     *cos(angleStartRotor+rad(90)+rad(90)))           calc(rMaxRotor     *sin(angleStartRotor+rad(90)+rad(90)))            zB) vlabel(BB3_090)

//Plane A:
//Bottom of rotor block#7 : sector without impeller 180-270 deg
(calc(rMaxRotor     *cos(angleStartRotor+rad(180)))                  calc(rMaxRotor     *sin(angleStartRotor+rad(180)))                   zA) vlabel(AA1_180)
(calc(rMaxRotor     *cos(angleStartRotor+rad(90)+rad(180)))          calc(rMaxRotor     *sin(angleStartRotor+rad(90)+rad(180)))           zA) vlabel(AA3_180)

//Plane B:
//Top of rotor block#7 : sector without impeller 180-270 deg
(calc(rMaxRotor     *cos(angleStartRotor+rad(180)))                  calc(rMaxRotor     *sin(angleStartRotor+rad(180)))                   zB) vlabel(BB1_180)
(calc(rMaxRotor     *cos(angleStartRotor+rad(90)+rad(180)))          calc(rMaxRotor     *sin(angleStartRotor+rad(90)+rad(180)))           zB) vlabel(BB3_180)

//Plane A:
//Bottom of rotor block#8 : sector without impeller 270-360 deg
(calc(rMaxRotor     *cos(angleStartRotor+rad(270)))                  calc(rMaxRotor     *sin(angleStartRotor+rad(270)))                   zA) vlabel(AA1_270)
(calc(rMaxRotor     *cos(angleStartRotor+rad(90)+rad(270)))          calc(rMaxRotor     *sin(angleStartRotor+rad(90)+rad(270)))           zA) vlabel(AA3_270)

//Plane B:
//Top of rotor block#8 : sector without impeller 270-360 deg
(calc(rMaxRotor     *cos(angleStartRotor+rad(270)))                  calc(rMaxRotor     *sin(angleStartRotor+rad(270)))                   zB) vlabel(BB1_270)
(calc(rMaxRotor     *cos(angleStartRotor+rad(90)+rad(270)))          calc(rMaxRotor     *sin(angleStartRotor+rad(90)+rad(270)))           zB) vlabel(BB3_270)

/////////////////////////////////
/////////////////////////////////

// Stator with impeller

//Plane C:
//Bottom of stator block#1 : sector without baffle 0-60 deg
(calc(rMinStator    *cos(angleStartStator))                           calc(rMinStator    *sin(angleStartStator))                          zC) vlabel(CC0_000)
(calc(rMinStator    *cos(angleStartStator+rad(60)))                   calc(rMinStator    *sin(angleStartStator+rad(60)))                  zC) vlabel(CC2_000)

//Plane D:
//Top of stator block#1 : sector without baffle 0-60 deg
(calc(rMinStator    *cos(angleStartStator))                           calc(rMinStator    *sin(angleStartStator))                          zD) vlabel(DD0_000)
(calc(rMinStator    *cos(angleStartStator+rad(60)))                   calc(rMinStator    *sin(angleStartStator+rad(60)))                  zD) vlabel(DD2_000)

//Plane C:
//Bottom of stator block#2 : sector without baffle 60-120 deg
(calc(rMinStator    *cos(angleStartStator+rad(60)))                   calc(rMinStator    *sin(angleStartStator+rad(60)))                  zC) vlabel(CC0_060)
(calc(rMinStator    *cos(angleStartStator+rad(60)+rad(60)))           calc(rMinStator    *sin(angleStartStator+rad(60)+rad(60)))          zC) vlabel(CC2_060)

//Plane D:
//Top of stator block#2 : sector without baffle 60-120 deg
(calc(rMinStator    *cos(angleStartStator+rad(60)))                   calc(rMinStator    *sin(angleStartStator+rad(60)))                  zD) vlabel(DD0_060)
(calc(rMinStator    *cos(angleStartStator+rad(60)+rad(60)))           calc(rMinStator    *sin(angleStartStator+rad(60)+rad(60)))          zD) vlabel(DD2_060)

//Plane C:
//Bottom of stator block#3 : sector without baffle 120-180 deg
(calc(rMinStator    *cos(angleStartStator+rad(120)))                  calc(rMinStator    *sin(angleStartStator+rad(120)))                 zC) vlabel(CC0_120)
(calc(rMinStator    *cos(angleStartStator+rad(60)+rad(120)))          calc(rMinStator    *sin(angleStartStator+rad(60)+rad(120)))         zC) vlabel(CC2_120)

//Plane D:
//Top of stator block#3 : sector without baffle 120-180 deg
(calc(rMinStator    *cos(angleStartStator+rad(120)))                  calc(rMinStator    *sin(angleStartStator+rad(120)))                 zD) vlabel(DD0_120)
(calc(rMinStator    *cos(angleStartStator+rad(60)+rad(120)))          calc(rMinStator    *sin(angleStartStator+rad(60)+rad(120)))         zD) vlabel(DD2_120)

//Plane C:
//Bottom of stator block#4 : sector without baffle 180-240 deg
(calc(rMinStator    *cos(angleStartStator+rad(180)))                  calc(rMinStator    *sin(angleStartStator+rad(180)))                 zC) vlabel(CC0_180)
(calc(rMinStator    *cos(angleStartStator+rad(60)+rad(180)))          calc(rMinStator    *sin(angleStartStator+rad(60)+rad(180)))         zC) vlabel(CC2_180)

//Plane D:
//Top of stator block#4 : sector without baffle 180-240 deg
(calc(rMinStator    *cos(angleStartStator+rad(180)))                  calc(rMinStator    *sin(angleStartStator+rad(180)))                 zD) vlabel(DD0_180)
(calc(rMinStator    *cos(angleStartStator+rad(60)+rad(180)))          calc(rMinStator    *sin(angleStartStator+rad(60)+rad(180)))         zD) vlabel(DD2_180)

//Plane C:
//Bottom of stator block#5 : sector without baffle 240-300 deg
(calc(rMinStator    *cos(angleStartStator+rad(240)))                  calc(rMinStator    *sin(angleStartStator+rad(240)))                 zC) vlabel(CC0_240)
(calc(rMinStator    *cos(angleStartStator+rad(60)+rad(240)))          calc(rMinStator    *sin(angleStartStator+rad(60)+rad(240)))         zC) vlabel(CC2_240)

//Plane D:
//Top of stator block#5 : sector without baffle 240-300 deg
(calc(rMinStator    *cos(angleStartStator+rad(240)))                  calc(rMinStator    *sin(angleStartStator+rad(240)))                 zD) vlabel(DD0_240)
(calc(rMinStator    *cos(angleStartStator+rad(60)+rad(240)))          calc(rMinStator    *sin(angleStartStator+rad(60)+rad(240)))         zD) vlabel(DD2_240)

//Plane C:
//Bottom of stator block#6 : sector without baffle 300-360 deg
(calc(rMinStator    *cos(angleStartStator+rad(300)))                  calc(rMinStator    *sin(angleStartStator+rad(300)))                 zC) vlabel(CC0_300)
(calc(rMinStator    *cos(angleStartStator+rad(60)+rad(300)))          calc(rMinStator    *sin(angleStartStator+rad(60)+rad(300)))         zC) vlabel(CC2_300)

//Plane D:
//Top of stator block#6 : sector without baffle 300-360 deg
(calc(rMinStator    *cos(angleStartStator+rad(300)))                  calc(rMinStator    *sin(angleStartStator+rad(300)))                 zD) vlabel(DD0_300)
(calc(rMinStator    *cos(angleStartStator+rad(60)+rad(300)))          calc(rMinStator    *sin(angleStartStator+rad(60)+rad(300)))         zD) vlabel(DD2_300)

//Plane C:
//Bottom of stator block#7 : sector with baffle 0-60 deg
(calc(rStatorBaffle*cos(angleStartStator))                            calc(rStatorBaffle*sin(angleStartStator))                           zC) vlabel(C0_000)
(calc(rMaxStator   *cos(angleStartStator+bafThick))                   calc(rMaxStator   *sin(angleStartStator+bafThick))                  zC) vlabel(C1_000)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)))                    calc(rStatorBaffle*sin(angleStartStator+rad(60)))                   zC) vlabel(C2_000)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks))          calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick))          zC) vlabel(C3_000)

//Plane D:
//Top of stator block#7 : sector with baffle 0-60 deg
(calc(rStatorBaffle*cos(angleStartStator))                            calc(rStatorBaffle*sin(angleStartStator))                           zD) vlabel(D0_000)
(calc(rMaxStator   *cos(angleStartStator+bafThick))                   calc(rMaxStator   *sin(angleStartStator+bafThick))                  zD) vlabel(D1_000)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)))                    calc(rStatorBaffle*sin(angleStartStator+rad(60)))                   zD) vlabel(D2_000)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks))          calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick))          zD) vlabel(D3_000)

//Plane C:
//Bottom of stator block#8 : sector with baffle 60-120 deg
(calc(rStatorBaffle*cos(angleStartStator+rad(60)))                    calc(rStatorBaffle*sin(angleStartStator+rad(60)))                   zC) vlabel(C0_060)
(calc(rMaxStator   *cos(angleStartStator+bafThick+rad(60)))           calc(rMaxStator   *sin(angleStartStator+bafThick+rad(60)))          zC) vlabel(C1_060)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)+rad(60)))            calc(rStatorBaffle*sin(angleStartStator+rad(60)+rad(60)))           zC) vlabel(C2_060)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks+rad(60)))  calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick+rad(60)))  zC) vlabel(C3_060)

//Plane D:
//Top of stator block#8 : sector with baffle 60-120 deg
(calc(rStatorBaffle*cos(angleStartStator+rad(60)))                    calc(rStatorBaffle*sin(angleStartStator+rad(60)))                   zD) vlabel(D0_060)
(calc(rMaxStator   *cos(angleStartStator+bafThick+rad(60)))           calc(rMaxStator   *sin(angleStartStator+bafThick+rad(60)))          zD) vlabel(D1_060)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)+rad(60)))            calc(rStatorBaffle*sin(angleStartStator+rad(60)+rad(60)))           zD) vlabel(D2_060)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks+rad(60)))  calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick+rad(60)))  zD) vlabel(D3_060)

//Plane C:
//Bottom of stator block#9 : sector with baffle 120-180 deg
(calc(rStatorBaffle*cos(angleStartStator+rad(120)))                   calc(rStatorBaffle*sin(angleStartStator+rad(120)))                  zC) vlabel(C0_120)
(calc(rMaxStator   *cos(angleStartStator+bafThick+rad(120)))          calc(rMaxStator   *sin(angleStartStator+bafThick+rad(120)))         zC) vlabel(C1_120)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)+rad(120)))           calc(rStatorBaffle*sin(angleStartStator+rad(60)+rad(120)))          zC) vlabel(C2_120)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks+rad(120))) calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick+rad(120))) zC) vlabel(C3_120)

//Plane D:
//Top of stator block#9 : sector with baffle 120-180 deg
(calc(rStatorBaffle*cos(angleStartStator+rad(120)))                   calc(rStatorBaffle*sin(angleStartStator+rad(120)))                  zD) vlabel(D0_120)
(calc(rMaxStator   *cos(angleStartStator+bafThick+rad(120)))          calc(rMaxStator   *sin(angleStartStator+bafThick+rad(120)))         zD) vlabel(D1_120)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)+rad(120)))           calc(rStatorBaffle*sin(angleStartStator+rad(60)+rad(120)))          zD) vlabel(D2_120)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks+rad(120))) calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick+rad(120))) zD) vlabel(D3_120)

//Plane C:
//Bottom of stator block#10 : sector with baffle 180-240 deg
(calc(rStatorBaffle*cos(angleStartStator+rad(180)))                   calc(rStatorBaffle*sin(angleStartStator+rad(180)))                  zC) vlabel(C0_180)
(calc(rMaxStator   *cos(angleStartStator+bafThick+rad(180)))          calc(rMaxStator   *sin(angleStartStator+bafThick+rad(180)))         zC) vlabel(C1_180)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)+rad(180)))           calc(rStatorBaffle*sin(angleStartStator+rad(60)+rad(180)))          zC) vlabel(C2_180)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks+rad(180))) calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick+rad(180))) zC) vlabel(C3_180)

//Plane D:
//Top of stator block#10 : sector with baffle 180-240 deg
(calc(rStatorBaffle*cos(angleStartStator+rad(180)))                   calc(rStatorBaffle*sin(angleStartStator+rad(180)))                  zD) vlabel(D0_180)
(calc(rMaxStator   *cos(angleStartStator+bafThick+rad(180)))          calc(rMaxStator   *sin(angleStartStator+bafThick+rad(180)))         zD) vlabel(D1_180)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)+rad(180)))           calc(rStatorBaffle*sin(angleStartStator+rad(60)+rad(180)))          zD) vlabel(D2_180)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks+rad(180))) calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick+rad(180))) zD) vlabel(D3_180)

//Plane C:
//Bottom of stator block#11 : sector with baffle 240-300 deg
(calc(rStatorBaffle*cos(angleStartStator+rad(240)))                   calc(rStatorBaffle*sin(angleStartStator+rad(240)))                  zC) vlabel(C0_240)
(calc(rMaxStator   *cos(angleStartStator+bafThick+rad(240)))          calc(rMaxStator   *sin(angleStartStator+bafThick+rad(240)))         zC) vlabel(C1_240)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)+rad(240)))           calc(rStatorBaffle*sin(angleStartStator+rad(60)+rad(240)))          zC) vlabel(C2_240)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks+rad(240))) calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick+rad(240))) zC) vlabel(C3_240)

//Plane D:
//Top of stator block#11 : sector with baffle 240-300 deg
(calc(rStatorBaffle*cos(angleStartStator+rad(240)))                   calc(rStatorBaffle*sin(angleStartStator+rad(240)))                  zD) vlabel(D0_240)
(calc(rMaxStator   *cos(angleStartStator+bafThick+rad(240)))          calc(rMaxStator   *sin(angleStartStator+bafThick+rad(240)))         zD) vlabel(D1_240)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)+rad(240)))           calc(rStatorBaffle*sin(angleStartStator+rad(60)+rad(240)))          zD) vlabel(D2_240)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks+rad(240))) calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick+rad(240))) zD) vlabel(D3_240)

//Plane C:
//Bottom of stator block#12 : sector with baffle 300-360 deg
(calc(rStatorBaffle*cos(angleStartStator+rad(300)))                   calc(rStatorBaffle*sin(angleStartStator+rad(300)))                  zC) vlabel(C0_300)
(calc(rMaxStator   *cos(angleStartStator+bafThick+rad(300)))          calc(rMaxStator   *sin(angleStartStator+bafThick+rad(300)))         zC) vlabel(C1_300)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)+rad(300)))           calc(rStatorBaffle*sin(angleStartStator+rad(60)+rad(300)))          zC) vlabel(C2_300)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks+rad(300))) calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick+rad(300))) zC) vlabel(C3_300)

//Plane D:
//Top of stator block#12 : sector with baffle 300-360 deg
(calc(rStatorBaffle*cos(angleStartStator+rad(300)))                   calc(rStatorBaffle*sin(angleStartStator+rad(300)))                  zD) vlabel(D0_300)
(calc(rMaxStator   *cos(angleStartStator+bafThick+rad(300)))          calc(rMaxStator   *sin(angleStartStator+bafThick+rad(300)))         zD) vlabel(D1_300)
(calc(rStatorBaffle*cos(angleStartStator+rad(60)+rad(300)))           calc(rStatorBaffle*sin(angleStartStator+rad(60)+rad(300)))          zD) vlabel(D2_300)
(calc(rMaxStator   *cos(angleStartStator+rad(60)-bafThicks+rad(300))) calc(rMaxStator   *sin(angleStartStator+rad(60)-bafThick+rad(300))) zD) vlabel(D3_300)

);

blocks
(
    // Rotor
    hex ( A0_000  A1_000   A1_090   A2_000  B0_000  B1_000   B1_090   B2_000 ) ( BLOCKSIZE_ROTOR_SECTOR_090 )  simpleGrading (1 1 grading)
    hex ( A0_090  A1_090   A1_180   A2_090  B0_090  B1_090   B1_180   B2_090 ) ( BLOCKSIZE_ROTOR_SECTOR_090 )  simpleGrading (1 1 grading)
    hex ( A0_180  A1_180   A1_270   A2_180  B0_180  B1_180   B1_270   B2_180 ) ( BLOCKSIZE_ROTOR_SECTOR_090 )  simpleGrading (1 1 grading)
    hex ( A0_270  A1_270   A1_000   A2_270  B0_270  B1_270   B1_000   B2_270 ) ( BLOCKSIZE_ROTOR_SECTOR_090 )  simpleGrading (1 1 grading)

    hex ( A1_000  AA1_000  AA1_090  A1_090  B1_000  BB1_000  BB1_090  B1_090 ) ( BLOCKSIZE_ROTOR_SECTOR_090 )  simpleGrading (1 1 grading)
    hex ( A1_090  AA1_090  AA1_180  A1_180  B1_090  BB1_090  BB1_180  B1_180 ) ( BLOCKSIZE_ROTOR_SECTOR_090 )  simpleGrading (1 1 grading)
    hex ( A1_180  AA1_180  AA1_270  A1_270  B1_180  BB1_180  BB1_270  B1_270 ) ( BLOCKSIZE_ROTOR_SECTOR_090 )  simpleGrading (1 1 grading)
    hex ( A1_270  AA1_270  AA1_000  A1_000  B1_270  BB1_270  BB1_000  B1_000 ) ( BLOCKSIZE_ROTOR_SECTOR_090 )  simpleGrading (1 1 grading)

    // Stator
    hex ( CC0_000 C0_000  C0_060  CC0_060 DD0_000 D0_000  D0_060  DD0_060 ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)
    hex ( CC0_060 C0_060  C0_120  CC0_120 DD0_060 D0_060  D0_120  DD0_120 ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)
    hex ( CC0_120 C0_120  C0_180  CC0_180 DD0_120 D0_120  D0_180  DD0_180 ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)
    hex ( CC0_180 C0_180  C0_240  CC0_240 DD0_180 D0_180  D0_240  DD0_240 ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)
    hex ( CC0_240 C0_240  C0_300  CC0_300 DD0_240 D0_240  D0_300  DD0_300 ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)
    hex ( CC0_300 C0_300  C0_000  CC0_000 DD0_300 D0_300  D0_000  DD0_000 ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)

    hex ( C0_000  C1_000  C3_000  C0_060  D0_000  D1_000  D3_000  D0_060  ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)
    hex ( C0_060  C1_060  C3_060  C0_120  D0_060  D1_060  D3_060  D0_120  ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)
    hex ( C0_120  C1_120  C3_120  C0_180  D0_120  D1_120  D3_120  D0_180  ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)
    hex ( C0_180  C1_180  C3_180  C0_240  D0_180  D1_180  D3_180  D0_240  ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)
    hex ( C0_240  C1_240  C3_240  C0_300  D0_240  D1_240  D3_240  D0_300  ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)
    hex ( C0_300  C1_300  C3_300  C0_000  D0_300  D1_300  D3_300  D0_000  ) ( BLOCKSIZE_STATOR_SECTOR_060 ) simpleGrading (1 1 grading)
);

edges
(
    // --- PLANE A: Bottom of rotor block#1
    arc  A0_000 A2_000   (calc(rMinRotor     *cos((angleStartRotor+rad(90))/2))          calc(rMinRotor     *sin((angleStartRotor+rad(90))/2))          zA)
    arc  A1_000 A1_090   (calc(rRotorImpeller*cos((angleStartRotor+rad(90))/2))          calc(rRotorImpeller*sin((angleStartRotor+rad(90))/2))          zA)

    // --- PLANE B: Top of rotor block#1
    arc  B0_000 B2_000   (calc(rMinRotor     *cos((angleStartRotor+rad(90))/2))          calc(rMinRotor     *sin((angleStartRotor+rad(90))/2))          zB)
    arc  B1_000 B1_090   (calc(rRotorImpeller*cos((angleStartRotor+rad(90))/2))          calc(rRotorImpeller*sin((angleStartRotor+rad(90))/2))          zB)

    // --- PLANE A: Bottom of rotor block#2
    arc  A0_090 A2_090   (calc(rMinRotor     *cos((angleStartRotor+rad(90))/2+rad(90)))  calc(rMinRotor     *sin((angleStartRotor+rad(90))/2+rad(90)))  zA)
    arc  A1_090 A1_180   (calc(rRotorImpeller*cos((angleStartRotor+rad(90))/2+rad(90)))  calc(rRotorImpeller*sin((angleStartRotor+rad(90))/2+rad(90)))  zA)

    // --- PLANE B: Top of rotor block#2
    arc  B0_090 B2_090   (calc(rMinRotor     *cos((angleStartRotor+rad(90))/2+rad(90)))  calc(rMinRotor     *sin((angleStartRotor+rad(90))/2+rad(90)))  zB)
    arc  B1_090 B1_180   (calc(rRotorImpeller*cos((angleStartRotor+rad(90))/2+rad(90)))  calc(rRotorImpeller*sin((angleStartRotor+rad(90))/2+rad(90)))  zB)

    // --- PLANE A: Bottom of rotor block#3
    arc  A0_180 A2_180   (calc(rMinRotor     *cos((angleStartRotor+rad(90))/2+rad(180))) calc(rMinRotor     *sin((angleStartRotor+rad(90))/2+rad(180))) zA)
    arc  A1_180 A1_270   (calc(rRotorImpeller*cos((angleStartRotor+rad(90))/2+rad(180))) calc(rRotorImpeller*sin((angleStartRotor+rad(90))/2+rad(180))) zA)

    // --- PLANE B: Top of rotor block#3
    arc  B0_180 B2_180   (calc(rMinRotor     *cos((angleStartRotor+rad(90))/2+rad(180))) calc(rMinRotor     *sin((angleStartRotor+rad(90))/2+rad(180))) zB)
    arc  B1_180 B1_270   (calc(rRotorImpeller*cos((angleStartRotor+rad(90))/2+rad(180))) calc(rRotorImpeller*sin((angleStartRotor+rad(90))/2+rad(180))) zB)

    // --- PLANE A: Bottom of rotor block#4
    arc  A0_270 A2_270   (calc(rMinRotor     *cos((angleStartRotor+rad(90))/2+rad(270))) calc(rMinRotor     *sin((angleStartRotor+rad(90))/2+rad(270))) zA)
    arc  A1_270 A1_000   (calc(rRotorImpeller*cos((angleStartRotor+rad(90))/2+rad(270))) calc(rRotorImpeller*sin((angleStartRotor+rad(90))/2+rad(270))) zA)

    // --- PLANE B: Top of rotor block#4
    arc  B0_270 B2_270   (calc(rMinRotor     *cos((angleStartRotor+rad(90))/2+rad(270))) calc(rMinRotor     *sin((angleStartRotor+rad(90))/2+rad(270))) zB)
    arc  B1_270 B1_000   (calc(rRotorImpeller*cos((angleStartRotor+rad(90))/2+rad(270))) calc(rRotorImpeller*sin((angleStartRotor+rad(90))/2+rad(270))) zB)

    // --- PLANE A: Bottom of rotor block#5
    arc  AA1_000 AA1_090 (calc(rMaxRotor     *cos((angleStartRotor+rad(90))/2))          calc(rMaxRotor     *sin((angleStartRotor+rad(90))/2))          zA)

    // --- PLANE B: Top of rotor block#5
    arc  BB1_000 BB1_090 (calc(rMaxRotor     *cos((angleStartRotor+rad(90))/2))          calc(rMaxRotor     *sin((angleStartRotor+rad(90))/2))          zB)

    // --- PLANE A: Bottom of rotor block#6
    arc  AA1_090 AA1_180 (calc(rMaxRotor     *cos((angleStartRotor+rad(90))/2+rad(90)))  calc(rMaxRotor     *sin((angleStartRotor+rad(90))/2+rad(90)))  zA)

    // --- PLANE B: Top of rotor block#6
    arc  BB1_090 BB1_180 (calc(rMaxRotor     *cos((angleStartRotor+rad(90))/2+rad(90)))  calc(rMaxRotor     *sin((angleStartRotor+rad(90))/2+rad(90)))  zB)

    // --- PLANE A: Bottom of rotor block#7
    arc  AA1_180 AA1_270 (calc(rMaxRotor     *cos((angleStartRotor+rad(90))/2+rad(180))) calc(rMaxRotor     *sin((angleStartRotor+rad(90))/2+rad(180))) zA)

    // --- PLANE B: Top of rotor block#7
    arc  BB1_180 BB1_270 (calc(rMaxRotor     *cos((angleStartRotor+rad(90))/2+rad(180))) calc(rMaxRotor     *sin((angleStartRotor+rad(90))/2+rad(180))) zB)

    // --- PLANE A: Bottom of rotor block#8
    arc  AA1_270 AA1_000 (calc(rMaxRotor     *cos((angleStartRotor+rad(90))/2+rad(270))) calc(rMaxRotor     *sin((angleStartRotor+rad(90))/2+rad(270))) zA)

    // --- PLANE B: Top of rotor block#8
    arc  BB1_270 BB1_000 (calc(rMaxRotor     *cos((angleStartRotor+rad(90))/2+rad(270))) calc(rMaxRotor     *sin((angleStartRotor+rad(90))/2+rad(270))) zB)

    // --- PLANE C: Bottom of stator block#1
    arc  CC0_000 CC0_060  (calc(rMinStator *cos((angleStartStator+rad(60))/2))           calc(rMinStator    *sin((angleStartStator+rad(60))/2))         zC)

    // --- PLANE D: Top of stator block#1
    arc  DD0_000 DD0_060  (calc(rMinStator *cos((angleStartStator+rad(60))/2))           calc(rMinStator    *sin((angleStartStator+rad(60))/2))         zD)

    // --- PLANE C: Bottom of stator block#2
    arc  CC0_060 CC0_120  (calc(rMinStator *cos((angleStartStator+rad(60))/2+rad(60)))   calc(rMinStator    *sin((angleStartStator+rad(60))/2+rad(60))) zC)

    // --- PLANE D: Top of stator block#2
    arc  DD0_060 DD0_120  (calc(rMinStator *cos((angleStartStator+rad(60))/2+rad(60)))   calc(rMinStator    *sin((angleStartStator+rad(60))/2+rad(60))) zD)

    // --- PLANE C: Bottom of stator block#3
    arc  CC0_120 CC0_180  (calc(rMinStator *cos((angleStartStator+rad(60))/2+rad(120)))  calc(rMinStator    *sin((angleStartStator+rad(60))/2+rad(120))) zC)

    // --- PLANE D: Top of stator block#3
    arc  DD0_120 DD0_180  (calc(rMinStator *cos((angleStartStator+rad(60))/2+rad(120)))  calc(rMinStator    *sin((angleStartStator+rad(60))/2+rad(120))) zD)

    // --- PLANE C: Bottom of stator block#4
    arc  CC0_180 CC0_240  (calc(rMinStator *cos((angleStartStator+rad(60))/2+rad(180)))  calc(rMinStator    *sin((angleStartStator+rad(60))/2+rad(180))) zC)

    // --- PLANE D: Top of stator block#4
    arc  DD0_180 DD0_240  (calc(rMinStator *cos((angleStartStator+rad(60))/2+rad(180)))  calc(rMinStator    *sin((angleStartStator+rad(60))/2+rad(180))) zD)

    // --- PLANE C: Bottom of stator block#5
    arc  CC0_240 CC0_300  (calc(rMinStator *cos((angleStartStator+rad(60))/2+rad(240)))  calc(rMinStator    *sin((angleStartStator+rad(60))/2+rad(240))) zC)

    // --- PLANE D: Top of stator block#5
    arc  DD0_240 DD0_300  (calc(rMinStator *cos((angleStartStator+rad(60))/2+rad(240)))  calc(rMinStator    *sin((angleStartStator+rad(60))/2+rad(240))) zD)

    // --- PLANE C: Bottom of stator block#6
    arc  CC0_300 CC0_000  (calc(rMinStator *cos((angleStartStator+rad(60))/2+rad(300)))  calc(rMinStator    *sin((angleStartStator+rad(60))/2+rad(300))) zC)

    // --- PLANE D: Top of stator block#6
    arc  DD0_300 DD0_000  (calc(rMinStator *cos((angleStartStator+rad(60))/2+rad(300)))  calc(rMinStator    *sin((angleStartStator+rad(60))/2+rad(300))) zD)

    // --- PLANE C: Bottom of stator block#7
    arc  C0_000 C0_060   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2))          calc(rStatorBaffle *sin((angleStartStator+rad(60))/2))          zC)
    arc  C1_000 C3_000   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2))          calc(rMaxStator    *sin((angleStartStator+rad(60))/2))          zC)

    // --- PLANE D: Top of stator block#7
    arc  D0_000 D0_060   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2))          calc(rStatorBaffle *sin((angleStartStator+rad(60))/2))          zD)
    arc  D1_000 D3_000   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2))          calc(rMaxStator    *sin((angleStartStator+rad(60))/2))          zD)

    // --- PLANE C: Bottom of stator block#8
    arc  C0_060 C0_120   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2+rad(60)))  calc(rStatorBaffle *sin((angleStartStator+rad(60))/2+rad(60)))  zC)
    arc  C1_060 C3_060   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2+rad(60)))  calc(rMaxStator    *sin((angleStartStator+rad(60))/2+rad(60)))  zC)

    // --- PLANE D: Top of stator block#8
    arc  D0_060 D0_120   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2+rad(60)))  calc(rStatorBaffle *sin((angleStartStator+rad(60))/2+rad(60)))  zD)
    arc  D1_060 D3_060   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2+rad(60)))  calc(rMaxStator    *sin((angleStartStator+rad(60))/2+rad(60)))  zD)

    // --- PLANE C: Bottom of stator block#9
    arc  C0_120 C0_180   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2+rad(120))) calc(rStatorBaffle *sin((angleStartStator+rad(60))/2+rad(120))) zC)
    arc  C1_120 C3_120   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2+rad(120))) calc(rMaxStator    *sin((angleStartStator+rad(60))/2+rad(120))) zC)

    // --- PLANE D: Top of stator block#9
    arc  D0_120 D0_180   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2+rad(180))) calc(rStatorBaffle *sin((angleStartStator+rad(60))/2+rad(180))) zD)
    arc  D1_120 D3_120   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2+rad(180))) calc(rMaxStator    *sin((angleStartStator+rad(60))/2+rad(180))) zD)

    // --- PLANE C: Bottom of stator block#10
    arc  C0_180 C0_240   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2+rad(240))) calc(rStatorBaffle *sin((angleStartStator+rad(60))/2+rad(240))) zC)
    arc  C1_180 C3_180   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2+rad(240))) calc(rMaxStator    *sin((angleStartStator+rad(60))/2+rad(240))) zC)

    // --- PLANE D: Top of stator block#10
    arc  D0_180 D0_240   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2+rad(240))) calc(rStatorBaffle *sin((angleStartStator+rad(60))/2+rad(240))) zD)
    arc  D1_180 D3_180   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2+rad(240))) calc(rMaxStator    *sin((angleStartStator+rad(60))/2+rad(240))) zD)

    // --- PLANE C: Bottom of stator block#11
    arc  C0_240 C0_300   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2+rad(300))) calc(rStatorBaffle *sin((angleStartStator+rad(60))/2+rad(300))) zC)
    arc  C1_240 C3_240   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2+rad(300))) calc(rMaxStator    *sin((angleStartStator+rad(60))/2+rad(300))) zC)

    // --- PLANE D: Top of stator block#11
    arc  D0_240 D0_300   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2+rad(300))) calc(rStatorBaffle *sin((angleStartStator+rad(60))/2+rad(300))) zD)
    arc  D1_240 D3_240   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2+rad(300))) calc(rMaxStator    *sin((angleStartStator+rad(60))/2+rad(300))) zD)

    // --- PLANE C: Bottom of stator block#12
    arc  C0_300 C0_000   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2+rad(360))) calc(rStatorBaffle *sin((angleStartStator+rad(60))/2+rad(360))) zC)
    arc  C1_300 C3_300   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2+rad(360))) calc(rMaxStator    *sin((angleStartStator+rad(60))/2+rad(360))) zC)

    // --- PLANE D: Top of stator block#12
    arc  D0_300 D0_000   (calc(rStatorBaffle*cos((angleStartStator+rad(60))/2+rad(360))) calc(rStatorBaffle *sin((angleStartStator+rad(60))/2+rad(360))) zD)
    arc  D1_300 D3_300   (calc(rMaxStator   *cos((angleStartStator+rad(60))/2+rad(360))) calc(rMaxStator    *sin((angleStartStator+rad(60))/2+rad(360))) zD)

);

boundary
(
    inlet
    {
        type patch;
        faces
        (
            ( A0_000 B0_000 B2_000 A2_000 )
            ( A0_090 B0_090 B2_090 A2_090 )
            ( A0_180 B0_180 B2_180 A2_180 )
            ( A0_270 B0_270 B2_270 A2_270 )
        );
    }

    outlet
    {
        type patch;
        faces
        (
            ( C1_000 C3_000 D3_000 D1_000 )
            ( C1_060 C3_060 D3_060 D1_060 )
            ( C1_120 C3_120 D3_120 D1_120 )
            ( C1_180 C3_180 D3_180 D1_180 )
            ( C1_240 C3_240 D3_240 D1_240 )
            ( C1_300 C3_300 D3_300 D1_300 )
        );
    }

    impellerWall
    {
        type wall;
        faces
        (
            ( A0_000 A1_000 B1_000 B0_000 )
            ( A0_090 A1_090 B1_090 B0_090 )
            ( A0_180 A1_180 B1_180 B0_180 )
            ( A0_270 A1_270 B1_270 B0_270 )

            ( A2_000 B2_000 B1_090 A1_090 )
            ( A2_090 B2_090 B1_180 A1_180 )
            ( A2_180 B2_180 B1_270 A1_270 )
            ( A2_270 B2_270 B1_000 A1_000 )

        );
    }

    baffleWall
    {
        type wall;
        faces
        (
            ( C0_000 C1_000 D1_000 D0_000 )
            ( C0_060 C1_060 D1_060 D0_060 )
            ( C0_120 C1_120 D1_120 D0_120 )
            ( C0_180 C1_180 D1_180 D0_180 )
            ( C0_240 C1_240 D1_240 D0_240 )
            ( C0_300 C1_300 D1_300 D0_300 )

            ( C0_060 D0_060 D3_000 C3_000 )
            ( C0_120 D0_120 D3_060 C3_060 )
            ( C0_180 D0_180 D3_120 C3_120 )
            ( C0_240 D0_240 D3_180 C3_180 )
            ( C0_300 D0_300 D3_240 C3_240 )
            ( C0_000 D0_000 D3_300 C3_300 )
        );
    }

    insideSlider
    {
        type ggi;
        shadowPatch     outsideSlider;
        zone            insideZone;
        bridgeOverlap   false;

        faces
        (
            ( AA1_000 AA1_090 BB1_090 BB1_000 )
            ( AA1_090 AA1_180 BB1_180 BB1_090 )
            ( AA1_180 AA1_270 BB1_270 BB1_180 )
            ( AA1_270 AA1_000 BB1_000 BB1_270 )
        );
    }

    outsideSlider
    {
        type            ggi;
        shadowPatch     insideSlider;
        zone            outsideZone;
        bridgeOverlap   false;

        faces
        (
            ( CC0_000 DD0_000 DD0_060 CC0_060 )
            ( CC0_060 DD0_060 DD0_120 CC0_120 )
            ( CC0_120 DD0_120 DD0_180 CC0_180 )
            ( CC0_180 DD0_180 DD0_240 CC0_240 )
            ( CC0_240 DD0_240 DD0_300 CC0_300 )
            ( CC0_300 DD0_300 DD0_000 CC0_000 )
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
