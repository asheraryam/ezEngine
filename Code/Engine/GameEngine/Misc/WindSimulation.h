#pragma once

#include <GameEngine/GameEngineDLL.h>

#include <Foundation/Containers/DynamicArray.h>
#include <Foundation/Time/Time.h>
#include <Foundation/Math/Declarations.h>

class EZ_GAMEENGINE_DLL ezWindSimulation
{
public:
  ezWindSimulation();
  ~ezWindSimulation();

  void Initialize(float fCellSize, ezUInt16 uiSizeX, ezUInt16 uiSizeY, ezUInt16 uiSizeZ = 1);
  void Step(ezTime tDelta);

  float GetCellSize() const { return m_fCellSize; }

  ezUInt16 GetSizeX() const { return m_uiSizeX; }
  ezUInt16 GetSizeY() const { return m_uiSizeY; }
  ezUInt16 GetSizeZ() const { return m_uiSizeZ; }

  EZ_ALWAYS_INLINE ezUInt32 Idx(ezUInt16 x, ezUInt16 y, ezUInt16 z = 0) const { return m_uiIndexOffsetZ * z + m_uiIndexOffsetY * y + x; }
  EZ_ALWAYS_INLINE bool IsVolumetric() const { return m_uiSizeZ > 1; }

  float* GetVelocitiesX() const { return m_pVelocities[0]; }
  float* GetVelocitiesY() const { return m_pVelocities[1]; }
  float* GetVelocitiesZ() const { return m_pVelocities[2]; }

  ezVec2 SampleVelocity2D(const ezVec2& vCellIdx) const;
  ezVec3 SampleVelocity3D(const ezVec3& vCellIdx) const;

private:
  void CopyPreviousVelocity(float* pDst, const float* pSrc);
  void LinearSolve(float* pDst, const float* pPrev);
  void Project2D(float* pDstU, float* pDstV, float* pScratch1, float* pScratch2);
  void Project3D(float* pDstU, float* pDstV, float* pDstW, float* pScratch1, float* pScratch2);
  void Advect(float* pDst, const float* pSrc);

  ezTime m_UpdateStep;
  float m_fCellSize = 1.0f;

  ezUInt16 m_uiSizeX = 0;
  ezUInt16 m_uiSizeY = 0;
  ezUInt16 m_uiSizeZ = 0;

  ezUInt32 m_uiNumCells = 0;

  ezUInt16 m_uiIndexOffsetY = 0;
  ezUInt16 m_uiIndexOffsetZ = 0;

  float* m_pVelocities[3] = {nullptr, nullptr, nullptr};
  float* m_pPrevVelocities[3] = {nullptr, nullptr, nullptr};

  ezDynamicArray<float> m_Values;
};
