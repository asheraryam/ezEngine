#include <PCH.h>
#include <EditorFramework/EditorApp/EditorApp.moc.h>
#include <Foundation/IO/OSFile.h>
#include <QStandardPaths>

ezString ezQtEditorApp::GetEditorDataFolder()
{
  ezStringBuilder sAppDir = ezOSFile::GetApplicationDirectory();
  sAppDir.AppendPath("../../../Data/Tools", ezUIServices::GetApplicationName());

  return sAppDir;
}

ezString ezQtEditorApp::GetPrecompiledToolsFolder()
{
  ezStringBuilder sPath = ezOSFile::GetApplicationDirectory();
  sPath.AppendPath("../../../Data/Tools/Precompiled");
  return sPath;

}

ezString ezQtEditorApp::GetDocumentDataFolder(const char* szDocument)
{
  ezStringBuilder sPath = szDocument;
  sPath.Append("_data");

  return sPath;
}


ezString ezQtEditorApp::GetEditorPreferencesFolder(bool bUserData)
{
  ezStringBuilder path = GetEditorDataFolder();
  path.AppendPath("Preferences");

  if (bUserData)
  {
    const QString userData = QStandardPaths::writableLocation(QStandardPaths::AppDataLocation);

    path = userData.toUtf8().data();
  }

  path.MakeCleanPath();
  return path;
}

ezString ezQtEditorApp::GetProjectPreferencesFolder(bool bUserData)
{
  ezStringBuilder path = ezToolsProject::GetSingleton()->GetProjectDataFolder();
  path.AppendPath("Preferences");

  if (bUserData)
  {
    const QString userData = QStandardPaths::writableLocation(QStandardPaths::AppDataLocation);

    ezStringBuilder ProjectName = ezToolsProject::GetSingleton()->GetProjectDirectory();
    ProjectName = ProjectName.GetFileName();

    path = userData.toUtf8().data();
    path.AppendPath("Projects", ProjectName);
  }

  path.MakeCleanPath();
  return path;
}

ezString ezQtEditorApp::GetDocumentPreferencesFolder(const ezDocument* pDocument, bool bUserData)
{
  ezStringBuilder path = GetDocumentDataFolder(pDocument->GetDocumentPath());
  path.AppendPath("Preferences");

  if (bUserData)
  {
    const ezStringBuilder sGuid = ezConversionUtils::ToString(pDocument->GetGuid());

    path = GetProjectPreferencesFolder(true);
    path.AppendPath(sGuid);
  }

  path.MakeCleanPath();
  return path;
}

