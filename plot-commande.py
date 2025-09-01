Remarque: j'ai installé la version gfortran 6.3 car la suivante 8.2 a un bug lors de la compilation avec l'impossibilité de compiler la fftw3 (problème de librairie)

rm out* restart-*;gfortran CMHD2D-v8.f90 -O3 -o CMHD2D-v8 -lfftw3;./CMHD2D-v8&

gfortran spectrum-anim.f90 -O2 -o spectrum-anim;./spectrum-anim
python3 plot-spectrum.py
gfortran spectrum-anim2.f90 -O2 -o spectrum-anim2;./spectrum-anim2
python3 plot-spectrum2.py
gfortran spectrum-anim3.f90 -O2 -o spectrum-anim3;./spectrum-anim3
python3 plot-spectrum3.py

