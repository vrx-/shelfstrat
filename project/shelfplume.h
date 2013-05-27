/*
** Options for shelf plume test case.
**
** Application flag:   SHELFPLUME
** Input script:       ocean_shelfplume.in
*/


#define ROMS_MODEL

#define UV_ADV
#define UV_COR
#define UV_LOGDRAG
#define TS_MPDATA

#define UV_PSOURCE
#define TS_PSOURCE

#define SALINITY
#define SOLVE3D
#define SPLINES
#define MASKING

#define DIAGNOSTICS_TS
#define DIAGNOSTICS_UV

#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS

#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_TGRADIENT
#define EAST_TNUDGING
#define EAST_M3GRADIENT

#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_TGRADIENT
#define WEST_M3GRADIENT

#define NORTH_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_TGRADIENT
#define NORTH_M3GRADIENT

#define GLS_MIXING
#define CANUTO_A
#define N2S2_HORAVG

#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_TOBC

#define ANA_INITIAL

#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX