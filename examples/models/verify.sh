#!/bin/bash

MODELS=()
# Loop through all .gms files in the current working directory
for file in *.gms; do
    # Check if the file actually exists (handles the case where no .gms files are present)
    [ -e "$file" ] || continue
    # Strip the '.gms' extension and append the base name to the array
    MODELS+=("${file%.gms}")
done
# Optional: Print the array to verify it worked
echo "Models found: ${MODELS[@]}"

# Tolerance for floating point comparison
TOLERANCE="0.0001"

echo "Starting GAMS objective value validation..."
echo "=========================================="

for model in "${MODELS[@]}"; do
    echo "Processing model: ${model}.gms"

    # 1. Run the model with the cuopt solver
    # We specify distinct output list files (o=...) so we can parse them
    echo "  -> Solving with cuopt..."
    gams "${model}.gms" \
         gdx="${model}_cuopt.gdx" \
         o="${model}_cuopt.lst" \
         lp=cuopt \
         mip=cuopt \
         qcp=cuopt \
         > "${model}_cuopt.log" 2>&1

    # 2. Run the model with the cplex solver
    echo "  -> Solving with cplex..."
    gams "${model}.gms" \
         gdx="${model}_cplex.gdx" \
         o="${model}_cplex.lst" \
         lp=cplex \
         mip=cplex \
         qcp=cplex \
         > "${model}_cplex.log" 2>&1

    # 3. Extract the Objective Values from the .lst files
    # We grab the last occurrence in case the model has multiple solves
    obj_cuopt=$(grep "OBJECTIVE VALUE" "${model}_cuopt.lst" | tail -n 1 | awk '{print $NF}')
    obj_cplex=$(grep "OBJECTIVE VALUE" "${model}_cplex.lst" | tail -n 1 | awk '{print $NF}')

    # Check if we successfully extracted the values (handles solve failures)
    if [ -z "$obj_cuopt" ] || [ -z "$obj_cplex" ]; then
        echo "  [!] WARNING: Could not extract objective value for ${model}.gms."
        echo "               cuopt obj: '${obj_cuopt}', cplex obj: '${obj_cplex}'"
        echo "               Check the .lst or .log files."
        echo "------------------------------------------"
        continue
    fi

    # 4. Compare the objective values using awk (Bash doesn't do floats natively)
    # Returns 1 if difference is greater than TOLERANCE, 0 otherwise
    is_diff=$(awk -v v1="$obj_cuopt" -v v2="$obj_cplex" -v tol="$TOLERANCE" '
        BEGIN {
            diff = v1 - v2;
            if (diff < 0) diff = -diff;
            if (diff > tol) print 1;
            else print 0;
        }')

    # 5. Check the comparison result and clean up on success
    if [ "$is_diff" -eq 1 ]; then
        echo "  [!] WARNING: Objectives differ!"
        echo "               cuopt: $obj_cuopt"
        echo "               cplex: $obj_cplex"
        echo "               Files retained for debugging."
    else
        echo "  [✓] SUCCESS: Objective values match (cuopt: $obj_cuopt, cplex: $obj_cplex)."
        echo "  -> Cleaning up generated files..."
        
        # Remove all generated files on success
        rm -f "${model}_cuopt.gdx" "${model}_cplex.gdx"
        rm -f "${model}_cuopt.lst" "${model}_cplex.lst"
        rm -f "${model}_cuopt.log" "${model}_cplex.log"
    fi
    echo "------------------------------------------"
done

echo "Validation complete."