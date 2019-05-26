#pragma once

#include <GameEngine/GameEngineDLL.h>

#include <Foundation/Containers/DynamicArray.h>
#include <Foundation/Math/Declarations.h>
#include <Foundation/Time/Time.h>

class EZ_GAMEENGINE_DLL ezWindSimulation
{
public:
  ezWindSimulation();
  ~ezWindSimulation();

  void Initialize(float fCellSize, ezUInt16 uiSizeX, ezUInt16 uiSizeY, ezUInt16 uiSizeZ = 1);
  void Step(ezTime tDelta);

  float GetCellSize() const { return m_fCellSize; }
  float GetInverseCellSize() const { return m_fInverseCellSize; }

  ezUInt16 GetSizeX() const { return m_uiSizeX; }
  ezUInt16 GetSizeY() const { return m_uiSizeY; }
  ezUInt16 GetSizeZ() const { return m_uiSizeZ; }

  EZ_ALWAYS_INLINE ezUInt32 Idx(ezUInt16 x, ezUInt16 y, ezUInt16 z = 0) const { return m_uiIndexOffsetZ * z + m_uiIndexOffsetY * y + x; }
  EZ_ALWAYS_INLINE bool IsVolumetric() const { return m_uiSizeZ > 1; }

  ezVec2* GetVelocities2D() const { return m_pVelocities2D[m_uiCurVelocities]; }
  ezVec3* GetVelocities3D() const { return m_pVelocities3D[m_uiCurVelocities]; }

  ezVec2 MapPositionToCellIdx(const ezVec2& vPosition) const;
  ezVec3 MapPositionToCellIdx(const ezVec3& vPosition) const;

  ezVec2 SampleVelocity2D(const ezVec2& vCellIdx) const;
  ezVec3 SampleVelocity3D(const ezVec3& vCellIdx) const;

private:
  ezTime m_UpdateStep;
  float m_fDampenFactor = 0.995f;
  float m_fCellSize = 1.0f;
  float m_fInverseCellSize = 1.0f;

  ezUInt16 m_uiSizeX = 0;
  ezUInt16 m_uiSizeY = 0;
  ezUInt16 m_uiSizeZ = 0;

  ezUInt16 m_uiIndexOffsetY = 0;
  ezUInt16 m_uiIndexOffsetZ = 0;

  ezUInt8 m_uiCurVelocities = 0;

  ezVec2* m_pVelocities2D[3] = {nullptr, nullptr, nullptr};
  ezVec2* m_pScratch2D = nullptr;
  ezVec3* m_pVelocities3D[3] = {nullptr, nullptr, nullptr};
  ezVec3* m_pScratch3D = nullptr;

  ezDynamicArray<float> m_Values;
};
