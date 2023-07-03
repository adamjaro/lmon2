
#set path to lmon2 executable in build

build="build"

export PATH=$PATH:`pwd`/$build:`pwd`/"reconstruction/python"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/$build:`pwd`/$build"/reconstruction"

