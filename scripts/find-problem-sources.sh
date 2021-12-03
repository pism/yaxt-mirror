#! /bin/bash
declare -a c_files f90_files
read -a c_files -d '' -r < <(find . ! -path ./config/checksrc/\* \
           ! -path ./config/workarounds/\* \
           ! -path ./include/cfortran.h \
           -type f \( -name \*.c -o -name \*.h \) -print)
for ((i=0; i < ${#c_files[@]}; ++i)) ; do
  f=${c_files[i]}
  declare -i problem_found
  problem_found=0
  if [[ "$f" != ./src/xt/*.h && "$f" != ./src/xt_quicksort_base.h ]] && ! grep '# *include *["<]config.h[>"]' "$f" >/dev/null ; then
    echo "$f does not contain include of config.h!" >&2
    problem_found=1
  fi
  if ! grep '\* *[@\\]file '"$(basename "$f")" "$f" >/dev/null ; then
    echo "$f contains no/malformed license/documentation header!" >&2
    problem_found=1
  fi
  if ! grep '^ *\* *Local Variables:' "$f" >/dev/null ; then
    echo "$f contains no/malformed local variable footer!" >&2
    problem_found=1
  fi
  if ((problem_found)); then
    ${VISUAL-${EDITOR-vi}} "$f"
    i=$((i-1))
  fi
done
read -a f90_files -d '' -r < <(find . ! -path ./config/checksrc/\* \
                       -type f \( -name \*.f90 \) -print)
for ((i=0; i < ${#f90_files[@]}; ++i)) ; do
  f=${f90_files[i]}
  declare -i problem_found
  problem_found=0
  if ! grep '#include *["<]fc_feature_defs.inc[>"]' "$f" >/dev/null ; then
    echo "$f does not contain include of fc_feature_defs.inc!" >&2
    problem_found=1
  fi
  if ! grep '![!>] *[@\\]file '"$(basename "$f")" "$f" >/dev/null ; then
    echo "$f contains no/malformed license/documentation header!" >&2
    problem_found=1
  fi
  if ! grep '^! Local Variables:' "$f" >/dev/null ; then
    echo "$f contains no/malformed local variable footer!" >&2
    problem_found=1
  fi
  if ((problem_found)); then
    ${VISUAL-${EDITOR-vi}} "$f"
    i=$((i-1))
  fi
done
