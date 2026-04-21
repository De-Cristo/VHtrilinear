#!/bin/bash

# Shared Stage 3 fixups for the VHtrilinear pipeline.
# This file is sourced by setup_env.sh and run_pipeline.sh.

stage3_require_file() {
    local file="$1"
    if [ ! -f "$file" ]; then
        echo "  [✗] Missing expected file: $file"
        return 1
    fi
}

stage3_patch_pdf_wrapper() {
    local file="$1"
    stage3_require_file "$file" || return 0

    if grep -q '^#include <boost/shared_ptr.hpp>' "$file"; then
        return 0
    fi

    sed -i '/^using namespace std;$/i #include <boost/shared_ptr.hpp>\n#include <boost/foreach.hpp>\n#include <boost/algorithm/string.hpp>' "$file"
}

stage3_patch_mg5_boost_wrappers() {
    local mg5dir="$1"
    local file

    while IFS= read -r -d '' file; do
        stage3_patch_pdf_wrapper "$file"
    done < <(find "$mg5dir" \( -path '*/Template/*/Source/PDF/pdf_lhapdf6.cc' -o -path '*/hz_*/Source/PDF/pdf_lhapdf6.cc' \) -type f -print0)
}

stage3_verify_container_boost_headers() {
    stage3_require_file /usr/include/boost/shared_ptr.hpp
    stage3_require_file /usr/include/boost/algorithm/string.hpp
}

stage3_patch_check_olp_pdflabel() {
    local file="$1"
    stage3_require_file "$file" || return 0
    sed -i "s/pdlabel='lhapdf'/pdlabel='nn23nlo'/g" "$file"
}

