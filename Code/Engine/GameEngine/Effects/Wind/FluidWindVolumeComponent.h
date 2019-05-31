#pragma once

#include <GameEngine/GameEngineDLL.h>

#include <Core/World/World.h>
#include <Core/World/Component.h>
#include <GameEngine/Effects/Wind/WindSimulation.h>
#include <GameEngine/Interfaces/WindWorldModule.h>

typedef ezComponentManagerSimple<class ezFluidWindVolumeComponent, ezComponentUpdateType::WhenSimulating> ezWindVolumeComponentManager;

class ezFluidWindVolume : public ezWindVolume
{
public:
  virtual ezResult GetWindAt(const ezVec3& vGlobalPosition, ezVec3& out_vWind) override;

  ezVec3 m_vPosition;
  ezQuat m_qRotation;
  ezWindSimulation m_Simulation;
};

class EZ_GAMEENGINE_DLL ezFluidWindVolumeComponent : public ezComponent
{
  EZ_DECLARE_COMPONENT_TYPE(ezFluidWindVolumeComponent, ezComponent, ezWindVolumeComponentManager);

public:
  ezFluidWindVolumeComponent();
  ~ezFluidWindVolumeComponent();

  void Update();

  virtual void SerializeComponent(ezWorldWriter& stream) const override;
  virtual void DeserializeComponent(ezWorldReader& stream) override;

  // ************************************* PROPERTIES ***********************************

  ezUInt8 m_uiCellsX = 32;
  ezUInt8 m_uiCellsY = 32;
  ezUInt8 m_uiCellsZ = 4;

  float m_fCellSize = 0.5f;
  bool m_bVisualize = false;

protected:
   
  virtual void OnSimulationStarted() override;
  virtual void OnDeactivated() override;

  void addStream(int x, int y, float force);
  void addDrop(int x, int y, float force);

  ezWindWorldModuleInterface* m_pWindModule = nullptr;
  ezFluidWindVolume m_FluidVolume;
};
