#include "GetFiles.h"

using namespace std;

const vector<string>& GetFiles::getFiles(const string& path)
{
    //文件句柄  
    intptr_t hFile = 0;
    //文件信息  
    struct _finddata_t fileinfo;
    string p;
    if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
    {
        do
        {
            //如果是目录,迭代之  
            //如果不是,加入列表  
            if ((fileinfo.attrib & _A_SUBDIR))
            {
                if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
                    getFiles(p.assign(path).append("\\").append(fileinfo.name));
            }
            else
            {
                m_FilesVector.push_back(p.assign(path).append("\\").append(fileinfo.name));
            }
        } while (_findnext(hFile, &fileinfo) == 0);
        _findclose(hFile);
    }
    return m_FilesVector;
}

#define MAX_PATH 80
#include <direct.h>
string GetFiles::GetWorkingDirectory()
{
    char buffer[MAX_PATH];
    _getcwd(buffer, MAX_PATH);
    return string(buffer);
}