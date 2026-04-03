Periodic sims.

- Run the 2D interpolation:
  ```
  python3 csv_to_interpolation2DTable_fixed_v2.py \
    --input slice.csv \
    --output UxProfile2D.table \
    --y-name Points:1 \
    --z-name Points:2 \
    --value-name U:0 \
    --dy 0.001 \
    --dz 0.001 \
    --mask-method circle \
    --circle-center-y 0 \
    --circle-center-z 0 \
    --circle-radius 0.05 \
    --fill-nearest
  ```
