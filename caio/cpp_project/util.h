#include <string>
#include <shlwapi.h>

#pragma once;

bool CopyDirTo(const std::wstring & source_folder, const std::wstring & target_folder )
{
    std::wstring new_sf = source_folder + L"\\*";
    WCHAR sf[MAX_PATH+1];
    WCHAR tf[MAX_PATH+1];

    wcscpy_s(sf, MAX_PATH, new_sf.c_str());
    wcscpy_s(tf, MAX_PATH, target_folder.c_str());

    sf[lstrlenW(sf)+1] = 0;
    tf[lstrlenW(tf)+1] = 0;

    SHFILEOPSTRUCTW s = { 0 };
    s.wFunc = FO_COPY;
    s.pTo = tf;
    s.pFrom = sf;
    s.fFlags = FOF_SILENT | FOF_NOCONFIRMMKDIR | FOF_NOCONFIRMATION | FOF_NOERRORUI | FOF_NO_UI;
    int res = SHFileOperationW( &s );

    return res == 0;
}