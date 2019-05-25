#include <GameEnginePCH.h>

#include <GameEngine/Interfaces/WindWorldModule.h>

// clang-format off
EZ_IMPLEMENT_WORLD_MODULE(ezWindWorldModuleInterface);

EZ_BEGIN_DYNAMIC_REFLECTED_TYPE(ezWindWorldModuleInterface, 1, ezRTTINoAllocator);
EZ_END_DYNAMIC_REFLECTED_TYPE;
// clang-format on

ezWindWorldModuleInterface::ezWindWorldModuleInterface(ezWorld* pWorld)
    : ezWorldModule(pWorld)
{
  m_vFallbackWind.SetZero();
}

 ezWindWorldModuleInterface::~ezWindWorldModuleInterface() = default;

ezVec3 ezWindWorldModuleInterface::GetWindAt(const ezVec3& vPosition)
{
  for (auto pVolume : m_WindVolumes)
  {
    ezVec3 result;
    if (pVolume->GetWindAt(vPosition, result).Succeeded())
      return result;
  }

  return m_vFallbackWind;
}

void ezWindWorldModuleInterface::SetFallbackWind(const ezVec3& vWind)
{
  m_vFallbackWind = vWind;
}

ezVec3 ezWindWorldModuleInterface::GetFallbackWind() const
{
  return m_vFallbackWind;
}

void ezWindWorldModuleInterface::AddWindVolume(ezWindVolume* pVolume)
{
  EZ_ASSERT_DEV(!m_WindVolumes.Contains(pVolume), "Cannot add the same wind volume twice");

  m_WindVolumes.PushBack(pVolume);
}

void ezWindWorldModuleInterface::RemoveWindVolume(ezWindVolume* pVolume)
{
  m_WindVolumes.RemoveAndSwap(pVolume);
}

EZ_STATICLINK_FILE(GameEngine, GameEngine_Interfaces_WindWorldModule);

