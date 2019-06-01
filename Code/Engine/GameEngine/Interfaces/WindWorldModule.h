#pragma once

#include <GameEngine/GameEngineDLL.h>

#include <Core/World/WorldModule.h>

class ezWindVolume
{
public:
  virtual ezResult GetWindAt(const ezVec3& vGlobalPosition, ezVec3& out_vWind) = 0;

  virtual void ApplyForceSphere(const ezVec3& vCenter, float fRadius, float fStrength) = 0;
};

class EZ_GAMEENGINE_DLL ezWindWorldModuleInterface : public ezWorldModule
{
  EZ_DECLARE_WORLD_MODULE();
  EZ_ADD_DYNAMIC_REFLECTION(ezWindWorldModuleInterface, ezWorldModule);

public:
  ezWindWorldModuleInterface(ezWorld* pWorld);
  ~ezWindWorldModuleInterface();

  virtual ezVec3 GetWindAt(const ezVec3& vGlobalPosition);

  virtual void SetFallbackWind(const ezVec3& vWind);
  virtual ezVec3 GetFallbackWind() const;

  void AddWindVolume(ezWindVolume* pVolume);
  void RemoveWindVolume(ezWindVolume* pVolume);

  virtual void ApplyForceSphere(const ezVec3& vCenter, float fRadius, float fStrength);

private:
  ezVec3 m_vFallbackWind;
  ezHybridArray<ezWindVolume*, 4> m_WindVolumes;
};
