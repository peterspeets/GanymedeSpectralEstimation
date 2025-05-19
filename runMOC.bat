set MOC_PATH=C:\Qt\6.9.0\mingw_64\bin\moc.exe
for %%f in (include\*.h) do (
    echo %%f
    "%MOC_PATH%" %%f -o moc_%%~nf.cpp
)