# Ux 
python3 csv_to_interpolation2DTable_fixed_v2.py --input slice.csv --output UxProfile2D.table --y-name Points:1 --z-name Points:2 --value-name U:0 --dy 0.001 --dz 0.001 --mask-method circle --circle-center-y 0 --circle-center-z 0 --circle-radius 0.05 --fill-nearest

# alpha.water
python3 csv_to_interpolation2DTable_fixed_v2.py --input slice.csv --output alphaWaterProfile2D.table --y-name Points:1 --z-name Points:2 --value-name alpha.water --dy 0.001 --dz 0.001 --mask-method circle --circle-center-y 0 --circle-center-z 0 --circle-radius 0.05 --fill-nearest

# k
python3 csv_to_interpolation2DTable_fixed_v2.py --input slice.csv --output kProfile2D.table --y-name Points:1 --z-name Points:2 --value-name k --dy 0.001 --dz 0.001 --mask-method circle --circle-center-y 0 --circle-center-z 0 --circle-radius 0.05 --fill-nearest

# nut
python3 csv_to_interpolation2DTable_fixed_v2.py --input slice.csv --output nutProfile2D.table --y-name Points:1 --z-name Points:2 --value-name nut --dy 0.001 --dz 0.001 --mask-method circle --circle-center-y 0 --circle-center-z 0 --circle-radius 0.05 --fill-nearest

# omega
python3 csv_to_interpolation2DTable_fixed_v2.py --input slice.csv --output omegaProfile2D.table --y-name Points:1 --z-name Points:2 --value-name omega --dy 0.001 --dz 0.001 --mask-method circle --circle-center-y 0 --circle-center-z 0 --circle-radius 0.05 --fill-nearest
