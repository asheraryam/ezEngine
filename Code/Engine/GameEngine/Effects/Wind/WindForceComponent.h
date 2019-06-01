#pragma once

#include <GameEngine/GameEngineDLL.h>

#include <Core/World/World.h>
#include <Core/World/Component.h>
#include <GameEngine/Interfaces/WindWorldModule.h>

class ezWindWorldModuleInterface;

typedef ezComponentManagerSimple<class ezWindForceComponent, ezComponentUpdateType::WhenSimulating> ezWindForceComponentManager;

class EZ_GAMEENGINE_DLL ezWindForceComponent : public ezComponent
{
  EZ_DECLARE_COMPONENT_TYPE(ezWindForceComponent, ezComponent, ezWindForceComponentManager);

public:
  ezWindForceComponent();
  ~ezWindForceComponent();

  void Update();

  virtual void SerializeComponent(ezWorldWriter& stream) const override;
  virtual void DeserializeComponent(ezWorldReader& stream) override;

  // ************************************* PROPERTIES ***********************************

  float m_fRadius = 1.0f;
  float m_fStrength = 1.0f;

protected:
  virtual void OnSimulationStarted() override;

  ezWindWorldModuleInterface* m_pWindModule = nullptr;
};
