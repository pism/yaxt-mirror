#! /bin/bash
declare -a c_files f90_files
c_files=()
f90_files=()

function add_file_to_category()
{
  case $1 in
    config/checksrc/*|config/work/arounds/*|include/cfortran.h)
      :
      ;;
    *.c|*.h)
      c_files+=("$1")
      ;;
    *.f90)
      f90_files+=("$1")
      ;;
  esac
}

if (( $# == 0 )); then
  read -a c_files -d '' -r < <(find . ! -path ./config/checksrc/\* \
           ! -path ./config/workarounds/\* \
           ! -path ./include/cfortran.h \
           -type f \( -name \*.c -o -name \*.h \) -print)
  read -a f90_files -d '' -r < <(find . ! -path ./config/checksrc/\* \
                       -type f \( -name \*.f90 \) -print)
elif [[ $# = 1 && x"$1" = x- ]]; then
  while read f ; do
    add_file_to_category "$f"
  done
else
  for f in "$@"; do
    add_file_to_category "$f"
  done
fi
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
