#include <PCH.h>

EZ_CREATE_SIMPLE_TEST_GROUP(Logging);

EZ_CREATE_SIMPLE_TEST(Logging, Log)
{
  ezLogWriter::HTML LogHTML;

  ezStringBuilder sOutputFolder1 = BUILDSYSTEM_OUTPUT_FOLDER;
  sOutputFolder1.AppendPath("FoundationTest");

  ezFileSystem::RegisterDataDirectoryFactory(ezDataDirectory::FolderType::Factory);
  ezFileSystem::AddDataDirectory(sOutputFolder1.GetData(), ezFileSystem::AllowWrites, "LoggingTest");

  LogHTML.BeginLog("Log_FoundationTest.htm", "FoundationTest");

  ezLog::AddLogWriter(ezLogWriter::Console::LogMessageHandler);
  ezLog::AddLogWriter(ezLogWriter::VisualStudio::LogMessageHandler);
  ezLog::AddLogWriter(ezLogWriter::HTML::LogMessageHandler, &LogHTML);

  ezStartup::PrintAllSubsystems();

  EZ_TEST_BLOCK(true, "Output")
  {
    EZ_LOG_BLOCK("Verse 1");

    ezLog::Success("%s", "This was a triumph.");
    ezLog::Info("%s", "I'm making a note here:");
    ezLog::Success("%s", "Huge Success");
    ezLog::Info("%s", "It's hard to overstate my satisfaction.");
    ezLog::Dev("%s", "Aperture Science. We do what we must, because we can,");
    ezLog::Debug("%s", "For the good of all of us, except the ones who are dead.");

    {
      EZ_LOG_BLOCK("Verse 2");

      ezLog::Dev("But there's no sense crying over every mistake.");
      ezLog::Debug("You just keep on trying 'till you run out of cake.");
      ezLog::Info("And the science gets done, and you make a neat gun");
      ezLog::DebugRegular("for the people who are still alive.");
    }

    {
      EZ_LOG_BLOCK("Verse 3");

      ezLog::Info("I'm not even angry.");
      ezLog::Debug("I'm being so sincere right now.");
      ezLog::Dev("Even though you broke my heart and killed me.");
      ezLog::Info("And tore me to pieces,");
      ezLog::Dev("and threw every piece into a fire.");
      ezLog::Info("As they burned it hurt because I was so happy for you.");
      ezLog::Dev("Now these points of data make a beautiful line");
      ezLog::Dev("and we're off the beta, we're releasing on time.");

      {
        EZ_LOG_BLOCK("Verse 4");

        ezLog::Info("So I'm glad I got burned,");
        ezLog::Debug("think of all the things we learned");
        ezLog::DebugRegular("for the people who are still alive.");
      
        {
          ezLogBlock b5("Verse 5");

          ezLog::Debug("Go ahead and leave me.");
          ezLog::Info("I think I prefer to stay inside.");
          ezLog::Dev("Maybe you'll find someone else, to help you.");
          ezLog::Dev("Maybe Black Mesa.");
          ezLog::Info("That was a joke. Haha. Fat chance.");
          ezLog::Success("Anyway, this cake is great.");
          ezLog::Success("It's so delicious and moist.");
          ezLog::Dev("Look at me still talking when there's science to do.");
          ezLog::Debug("When I look up there it makes me glad I'm not you.");
          ezLog::Info("I've experiments to run,");
          ezLog::DebugRegular("there is research to be done on the people who are still alive.");
        }
      }
    }
  }
    
  {
    EZ_LOG_BLOCK("Verse 6");

    ezLog::Dev("And believe me I am still alive.");
    ezLog::Info("I'm doing science and I'm still alive.");
    ezLog::Success("I feel fantastic and I'm still alive.");
    ezLog::Dev("While you're dying I'll be still alive.");
    ezLog::Dev("And when you're dead I will be, still alive.");
    ezLog::DebugRegular("Still alive, still alive.");
  }

  ezLog::RemoveLogWriter(ezLogWriter::Console::LogMessageHandler);
  ezLog::RemoveLogWriter(ezLogWriter::VisualStudio::LogMessageHandler);
  ezLog::RemoveLogWriter(ezLogWriter::HTML::LogMessageHandler, &LogHTML);

  LogHTML.EndLog();

  ezFileSystem::RemoveDataDirectoryGroup("LoggingTest");
  ezFileSystem::ClearAllDataDirectoryFactories();
}


