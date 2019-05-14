#pragma once

#include <GameEngine/GameEngineDLL.h>

#include <Core/World/World.h>
#include <Core/World/Component.h>
#include <GameEngine/Misc/WindSimulation.h>

typedef ezComponentManagerSimple<class ezWindVolumeComponent, ezComponentUpdateType::WhenSimulating> ezWindVolumeComponentManager;

class EZ_GAMEENGINE_DLL ezWindVolumeComponent : public ezComponent
{
  EZ_DECLARE_COMPONENT_TYPE(ezWindVolumeComponent, ezComponent, ezWindVolumeComponentManager);

public:
  ezWindVolumeComponent();

  void Update();

  virtual void SerializeComponent(ezWorldWriter& stream) const override;
  virtual void DeserializeComponent(ezWorldReader& stream) override;

  // ************************************* PROPERTIES ***********************************

protected:
   
  virtual void OnSimulationStarted() override;

  void addStream(int x, int y, float force);
  void addDrop(int x, int y, float force);

  ezWindSimulation m_Wind;
};
