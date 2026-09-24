#!/bin/sh
# This script generates the abinit documentation in ROBODOC format in the directory tmp-robodoc
echo "Will generate ROBODOC documentation in tmp-robodoc (requires robodoc)"

# Check the dependency before removing any previously generated documentation.
if ! command -v robodoc >/dev/null 2>&1; then
    echo "Error: robodoc is required but was not found in PATH." >&2
    echo "Install ROBODoc or add its executable directory to PATH, then run this script again." >&2
    exit 127
fi

#rm -rf tmp-robodoc robodoc-html && mkdir tmp-robodoc
#cp -rf ./src/[0-9]* tmp-robodoc
#cp ./config/robodoc/robodoc-html.rc tmp-robodoc/robodoc.rc
#cd tmp-robodoc && rm -f */*.in && rm -f */interfaces* && robodoc > ../robodoc.log 2> ../robodoc.err
#exit_status=`cat ../robodoc.err | wc -l`
#if test $exit_status -ne 0 ; then
#  cat ../doc/developers/robodoc.doc.txt >> robodoc.err
#  cat ../robodoc.err
#fi
#
#echo "Exit status: " $exit_status
#exit $exit_status

rm -rf tmp-robodoc robodoc-html && mkdir tmp-robodoc
# 39_libpaw is a symlink to shared/libpaw/src and is copied explicitly below.
# Excluding the link here prevents every LibPAW header from being staged twice.
for src_dir in ./shared/common/src/[0-3]*; do
    if test "$(basename "$src_dir")" != "39_libpaw"; then
        cp -rf "$src_dir" tmp-robodoc
    fi
done
cp -rf ./shared/libpaw/src tmp-robodoc/39_libpaw
cp -rf ./src/[4-9]* tmp-robodoc
cp ./config/robodoc/robodoc-html.rc tmp-robodoc/robodoc.rc
cd tmp-robodoc && rm -f */*.in && robodoc > ../robodoc.log 2> ../robodoc.err && cd ..
exit_status=`cat robodoc.err | wc -l`

# ROBODoc has no option to name its entry page "index.html" -- --multidoc
# --index always names it masterindex.html, and --toc always names the
# table-of-contents page toc_index.html (already the page abibuildbot's
# AbiInfo.robodoc_url links to from the build page, see mysteps.py). Without
# an index.html, a plain visit to the published site's root (e.g.
# https://robodoc.abinit.org/) hits no file nginx recognizes as a default
# index and falls back to a raw directory listing instead of any
# ROBODoc-generated page. Copy toc_index.html to index.html so the root URL
# lands on the same page the build-page link already uses, instead of
# picking a third, inconsistent entry point.
robodoc_doc_dir=tmp-robodoc/www/robodoc
if test -f "$robodoc_doc_dir/toc_index.html" ; then
    cp -f "$robodoc_doc_dir/toc_index.html" "$robodoc_doc_dir/index.html"
else
    echo "Error: $robodoc_doc_dir/toc_index.html was not generated -- cannot create index.html" >&2
    exit 1
fi
#mv -f tmp-robodoc/www/robodoc robodoc-html
#tardir=robodoc-html && tar --format=ustar -chf - "$tardir" | GZIP=--best gzip -c >robodoc-html-8.11.8.tar.gz
#rm -rf robodoc-html tmp-robodoc
#cat ./doc/developers/robodoc.doc.txt >> robodoc.err

if test $exit_status -ne 0 ; then
    cat robodoc.err
fi

echo "Exit status: " $exit_status
exit $exit_status
