#!/usr/bin/perl
$latex         = 'platex %O -src-specials -shell-escape -interaction=nonstopmode -synctex=1 -kanji=utf8 %S';
$bibtex        = 'pbibtex %O %B -kanji=utf8';
$dvipdf        = 'dvipdfmx %O -o %D %S';
$pdf_mode      = 3; # use dvipdfmx

# Use SumatraPDF and atom inverse search
$pdf_previewer = '"C:\Program Files\SumatraPDF\SumatraPDF.exe" -reuse-instance -inverse-search "\"%USERPROFILE%\AppData\Local\atom\bin\atom.cmd\" \"%f:%l\"" %O %S';

# if you do not need to preview in -pvc option (preview document and countinuously update mode)
# $pdf_previewer = 'exit';
