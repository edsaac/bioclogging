cd $1
foamExec $2
foamToVTK
cd postProcessing
$HOME/NUmodel/bioclogging/src/scripts/pointFiles.sh boundaryProbes Cw
cd ..
paraview
cd ..
