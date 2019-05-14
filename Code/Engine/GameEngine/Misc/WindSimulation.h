#pragma once

#include <GameEngine/GameEngineDLL.h>

#include <Foundation/Time/Time.h>
#include <Foundation/Containers/DynamicArray.h>

class EZ_GAMEENGINE_DLL ezWindSimulation
{
public:
  ezWindSimulation();
  ~ezWindSimulation();

  void Initialize(ezUInt16 uiSizeX, ezUInt16 uiSizeY, ezUInt16 uiSizeZ = 1);
  void Step();

  // private:
  EZ_ALWAYS_INLINE ezUInt32 Idx(ezUInt16 x, ezUInt16 y, ezUInt16 z = 0) const { return m_uiIndexOffsetZ * z + m_uiIndexOffsetY * y + x; }
  EZ_ALWAYS_INLINE bool IsVolumetric() const { return m_uiSizeZ > 1; }

  void AddTimeScaled(float* pDst, const float* pSrc);
  void ClearBounds(float* pDst);
  void LinearSolve(float* pDst, const float* pPrev);
  void Project(float* pDstU, float* pDstV, float* pScratchU, float* pScratchV);
  void Advect(float* pDst, const float* pSrc);

  ezTime m_UpdateStep;

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
