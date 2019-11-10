MATLAB="/Applications/MATLAB_R2017a.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/torbjri/Library/Application Support/MathWorks/MATLAB/R2017a"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for eigenValueSphere3" > eigenValueSphere3_mex.mki
echo "CC=$CC" >> eigenValueSphere3_mex.mki
echo "CFLAGS=$CFLAGS" >> eigenValueSphere3_mex.mki
echo "CLIBS=$CLIBS" >> eigenValueSphere3_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> eigenValueSphere3_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> eigenValueSphere3_mex.mki
echo "CXX=$CXX" >> eigenValueSphere3_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> eigenValueSphere3_mex.mki
echo "CXXLIBS=$CXXLIBS" >> eigenValueSphere3_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> eigenValueSphere3_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> eigenValueSphere3_mex.mki
echo "LDFLAGS=$LDFLAGS" >> eigenValueSphere3_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> eigenValueSphere3_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> eigenValueSphere3_mex.mki
echo "Arch=$Arch" >> eigenValueSphere3_mex.mki
echo "LD=$LD" >> eigenValueSphere3_mex.mki
echo OMPFLAGS= >> eigenValueSphere3_mex.mki
echo OMPLINKFLAGS= >> eigenValueSphere3_mex.mki
echo "EMC_COMPILER=clang" >> eigenValueSphere3_mex.mki
echo "EMC_CONFIG=optim" >> eigenValueSphere3_mex.mki
"/Applications/MATLAB_R2017a.app/bin/maci64/gmake" -B -f eigenValueSphere3_mex.mk
