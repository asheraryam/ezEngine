#include <GameEnginePCH.h>

#include <Foundation/Math/Vec3.h>
#include <Foundation/Profiling/Profiling.h>
#include <Foundation/Threading/DelegateTask.h>
#include <Foundation/Threading/TaskSystem.h>
#include <GameEngine/Misc/WindSimulation.h>

#if 1

// uses the empty boundary to affect the simulation

#  define ClampMinX 0
#  define ClampMaxX (uiSizeX + 1)
#  define ClampMinY 0
#  define ClampMaxY (uiSizeY + 1)
#  define ClampMinZ 0
#  define ClampMaxZ (uiSizeZ + 1)
#  define ClampOffset 0.0f

#else

// clamps sampling points to the inner values (no boundary)

#  define ClampMinX 1
#  define ClampMaxX (uiSizeX)
#  define ClampMinY 1
#  define ClampMaxY (uiSizeY)
#  define ClampMinZ 1
#  define ClampMaxZ (uiSizeZ)
#  define ClampOffset 1.0f

#endif

#define PrepareIdx2D const ezUInt32 uiIdxOffsetY = uiSizeX + 2;
#define PrepareIdx3D PrepareIdx2D const ezUInt32 uiIdxOffsetZ = (uiSizeY + 2) * (uiSizeX + 2);

#define Idx2D(x, y) (uiIdxOffsetY * y + x)
#define Idx3D(x, y, z) (uiIdxOffsetZ * z + uiIdxOffsetY * y + x)

namespace
{
  void CopyPreviousVelocity(ezVec2* pDst, const ezVec2* pSrc, ezVec2* pInputs, ezUInt32 uiNumCells, float fDampenFactor)
  {
    for (ezUInt32 i = 0; i < uiNumCells; ++i)
    {
      pDst[i] = pSrc[i] * fDampenFactor + pInputs[i];
    }

    ezMemoryUtils::ZeroFill(pInputs, uiNumCells);
  }

  void CopyPreviousVelocity(ezVec3* pDst, const ezVec3* pSrc, ezVec3* pInputs, ezUInt32 uiNumCells, float fDampenFactor)
  {
    for (ezUInt32 i = 0; i < uiNumCells; ++i)
    {
      pDst[i] = pSrc[i] * fDampenFactor + pInputs[i];
    }

    ezMemoryUtils::ZeroFill(pInputs, uiNumCells);
  }

  void LinearSolve2D(float* pDst, const float* pPrev, const ezUInt16 uiSizeX, const ezUInt16 uiSizeY)
  {
    const ezUInt32 uiNumIterations = 20;
    const float fNorm = 1.0f / 4.0f;

    PrepareIdx2D;

    for (ezUInt32 k = 0; k < uiNumIterations; k++)
    {
      for (ezUInt32 y = 1; y <= uiSizeY; ++y)
      {
        const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
        const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

        for (ezUInt32 x = 1; x <= uiSizeX; ++x)
        {
          const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
          const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

          pDst[Idx2D(x, y)] =
            (pPrev[Idx2D(x, y)] + (pDst[Idx2D(xm, y)] + pDst[Idx2D(xp, y)] + pDst[Idx2D(x, ym)] + pDst[Idx2D(x, yp)])) * fNorm;
        }
      }
    }
  }

  void LinearSolve3D(float* pDst, const float* pPrev, const ezUInt16 uiSizeX, const ezUInt16 uiSizeY, const ezUInt16 uiSizeZ)
  {
    const ezUInt32 uiNumIterations = 20;
    const float fNorm = 1.0f / 6.0f;

    PrepareIdx3D;

    for (ezUInt32 k = 0; k < uiNumIterations; k++)
    {
      for (ezInt32 z = 1; z <= uiSizeZ; ++z)
      {
        const ezInt32 zm = ezMath::Max<ezInt32>(z - 1, ClampMinZ);
        const ezInt32 zp = ezMath::Min<ezInt32>(z + 1, ClampMaxZ);

        for (ezInt32 y = 1; y <= uiSizeY; ++y)
        {
          const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
          const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

          for (ezInt32 x = 1; x <= uiSizeX; ++x)
          {
            const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
            const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

            pDst[Idx3D(x, y, z)] = (pPrev[Idx3D(x, y, z)] + (pDst[Idx3D(xm, y, z)] + pDst[Idx3D(xp, y, z)] + pDst[Idx3D(x, ym, z)] +
                                                              pDst[Idx3D(x, yp, z)] + pDst[Idx3D(x, y, zm)] + pDst[Idx3D(x, y, zp)])) *
                                   fNorm;
          }
        }
      }
    }
  }

  void Project2D(ezVec2* pSrcDst, float* pScratch1, float* pScratch2, const ezUInt16 uiSizeX, const ezUInt16 uiSizeY)
  {
    const float fNorm = -0.5f;
    const ezUInt32 uiNumCells = (uiSizeX + 2) * (uiSizeY + 2);

    PrepareIdx2D;

    for (ezUInt32 y = 1; y <= uiSizeY; ++y)
    {
      const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
      const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

      for (ezUInt32 x = 1; x <= uiSizeX; ++x)
      {
        const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
        const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

        pScratch2[Idx2D(x, y)] =
          (pSrcDst[Idx2D(xp, y)].x - pSrcDst[Idx2D(xm, y)].x + pSrcDst[Idx2D(x, yp)].y - pSrcDst[Idx2D(x, ym)].y) * fNorm;
      }
    }

    ezMemoryUtils::ZeroFill(pScratch1, uiNumCells);

    LinearSolve2D(pScratch1, pScratch2, uiSizeX, uiSizeY);

    const float fNorm2 = 0.5f;

    for (ezUInt32 y = 1; y <= uiSizeY; ++y)
    {
      const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
      const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

      for (ezUInt32 x = 1; x <= uiSizeX; ++x)
      {
        const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
        const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

        ezVec2 sub;
        sub.x = fNorm2 * (pScratch1[Idx2D(xp, y)] - pScratch1[Idx2D(xm, y)]);
        sub.y = fNorm2 * (pScratch1[Idx2D(x, yp)] - pScratch1[Idx2D(x, ym)]);

        pSrcDst[Idx2D(x, y)] -= sub;
      }
    }
  }

  void Project3D(
    ezVec3* pSrcDst, float* pScratch1, float* pScratch2, const ezUInt16 uiSizeX, const ezUInt16 uiSizeY, const ezUInt16 uiSizeZ)
  {
    const float fNorm = -0.33f;
    const ezUInt32 uiNumCells = (uiSizeX + 2) * (uiSizeY + 2) * (uiSizeZ + 2);

    PrepareIdx3D;

    for (ezUInt32 z = 1; z <= uiSizeZ; ++z)
    {
      const ezInt32 zm = ezMath::Max<ezInt32>(z - 1, ClampMinZ);
      const ezInt32 zp = ezMath::Min<ezInt32>(z + 1, ClampMaxZ);

      for (ezUInt32 y = 1; y <= uiSizeY; ++y)
      {
        const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
        const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

        for (ezUInt32 x = 1; x <= uiSizeX; ++x)
        {
          const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
          const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

          pScratch2[Idx3D(x, y, z)] = fNorm * (pSrcDst[Idx3D(xp, y, z)].x - pSrcDst[Idx3D(xm, y, z)].x + pSrcDst[Idx3D(x, yp, z)].y -
                                                pSrcDst[Idx3D(x, ym, z)].y + pSrcDst[Idx3D(x, y, zp)].z - pSrcDst[Idx3D(x, y, zm)].z);
        }
      }
    }

    ezMemoryUtils::ZeroFill(pScratch1, uiNumCells);

    LinearSolve3D(pScratch1, pScratch2, uiSizeX, uiSizeY, uiSizeZ);

    const float fNorm2 = 0.33f;

    for (ezUInt32 z = 1; z <= uiSizeZ; ++z)
    {
      const ezInt32 zm = ezMath::Max<ezInt32>(z - 1, ClampMinZ);
      const ezInt32 zp = ezMath::Min<ezInt32>(z + 1, ClampMaxZ);

      for (ezUInt32 y = 1; y <= uiSizeY; ++y)
      {
        const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
        const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

        for (ezUInt32 x = 1; x <= uiSizeX; ++x)
        {
          const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
          const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

          ezVec3 sub;
          sub.x = fNorm2 * (pScratch1[Idx3D(xp, y, z)] - pScratch1[Idx3D(xm, y, z)]);
          sub.y = fNorm2 * (pScratch1[Idx3D(x, yp, z)] - pScratch1[Idx3D(x, ym, z)]);
          sub.z = fNorm2 * (pScratch1[Idx3D(x, y, zp)] - pScratch1[Idx3D(x, y, zm)]);

          pSrcDst[Idx3D(x, y, z)] -= sub;
        }
      }
    }
  }

  void Advect2D(ezVec2* pDst, const ezVec2* pSrc, const ezUInt16 uiSizeX, const ezUInt16 uiSizeY, const float deltaTime)
  {
    const float minX = 0.5f + ClampOffset;
    const float minY = 0.5f + ClampOffset;
    const float maxX = uiSizeX + 0.5f - ClampOffset;
    const float maxY = uiSizeY + 0.5f - ClampOffset;

    PrepareIdx2D;

    for (ezUInt32 j = 1; j <= uiSizeY; ++j)
    {
      for (ezUInt32 i = 1; i <= uiSizeX; ++i)
      {
        // compute reverse velocity sample position
        const ezVec2 samplePos = pSrc[Idx2D(i, j)];
        const float x = ezMath::Clamp(i - deltaTime * samplePos.x, minX, maxX);
        const float y = ezMath::Clamp(j - deltaTime * samplePos.y, minY, maxY);

        // bilinear interpolate from the 4 sample position cells
        const ezUInt32 i0 = (ezUInt32)x;
        const ezUInt32 i1 = i0 + 1;
        const ezUInt32 j0 = (ezUInt32)y;
        const ezUInt32 j1 = j0 + 1;

        const float s1 = x - i0;
        const float s0 = 1 - s1;
        const float t1 = y - j0;
        const float t0 = 1 - t1;

        pDst[Idx2D(i, j)] =
          s0 * (t0 * pSrc[Idx2D(i0, j0)] + t1 * pSrc[Idx2D(i0, j1)]) + s1 * (t0 * pSrc[Idx2D(i1, j0)] + t1 * pSrc[Idx2D(i1, j1)]);
      }
    }
  }

  void Advect3D(
    ezVec3* pDst, const ezVec3* pSrc, const ezUInt16 uiSizeX, const ezUInt16 uiSizeY, const ezUInt16 uiSizeZ, const float deltaTime)
  {
    const float minX = 0.5f + ClampOffset;
    const float minY = 0.5f + ClampOffset;
    const float minZ = 0.5f + ClampOffset;
    const float maxX = uiSizeX + 0.5f - ClampOffset;
    const float maxY = uiSizeY + 0.5f - ClampOffset;
    const float maxZ = uiSizeZ + 0.5f - ClampOffset;

    PrepareIdx3D;

    for (ezUInt32 k = 1; k <= uiSizeZ; ++k)
    {
      for (ezUInt32 j = 1; j <= uiSizeY; ++j)
      {
        for (ezUInt32 i = 1; i <= uiSizeX; ++i)
        {
          // compute reverse velocity sample position
          const ezVec3 samplePos = pSrc[Idx3D(i, j, k)];
          const float x = ezMath::Clamp(i - deltaTime * samplePos.x, minX, maxX);
          const float y = ezMath::Clamp(j - deltaTime * samplePos.y, minY, maxY);
          const float z = ezMath::Clamp(k - deltaTime * samplePos.z, minZ, maxZ);

          // trilinear interpolate from the 8 sample position cells
          const ezUInt32 i0 = (ezUInt32)x;
          const ezUInt32 i1 = i0 + 1;
          const ezUInt32 j0 = (ezUInt32)y;
          const ezUInt32 j1 = j0 + 1;
          const ezUInt32 k0 = (ezUInt32)z;
          const ezUInt32 k1 = k0 + 1;

          const float s1 = x - i0;
          const float s0 = 1 - s1;
          const float t1 = y - j0;
          const float t0 = 1 - t1;
          const float r1 = z - k0;
          const float r0 = 1 - r1;

          pDst[Idx3D(i, j, k)] = r0 * (s0 * (t0 * pSrc[Idx3D(i0, j0, k0)] + t1 * pSrc[Idx3D(i0, j1, k0)]) +
                                        s1 * (t0 * pSrc[Idx3D(i1, j0, k0)] + t1 * pSrc[Idx3D(i1, j1, k0)])) +
                                 r1 * (s0 * (t0 * pSrc[Idx3D(i0, j0, k1)] + t1 * pSrc[Idx3D(i0, j1, k1)]) +
                                        s1 * (t0 * pSrc[Idx3D(i1, j0, k1)] + t1 * pSrc[Idx3D(i1, j1, k1)]));
        }
      }
    }
  }

  void StepWindSimulation2D(float deltaTime, float fDampenFactor, ezUInt16 uiSizeX, ezUInt16 uiSizeY, const ezVec2* pSrc, ezVec2* pDst,
    ezVec2* pInputs, ezVec2* pScratch)
  {
    const ezUInt32 uiNumCells = (uiSizeX + 2) * (uiSizeY + 2);

    CopyPreviousVelocity(pScratch, pSrc, pInputs, uiNumCells, fDampenFactor);

    Project2D(pScratch, pDst->GetData(), pDst->GetData() + uiNumCells, uiSizeX, uiSizeY);
    Advect2D(pDst, pScratch, uiSizeX, uiSizeY, deltaTime);
    Project2D(pDst, pScratch->GetData(), pScratch->GetData() + uiNumCells, uiSizeX, uiSizeY);
  }

  void StepWindSimulation3D(float deltaTime, float fDampenFactor, ezUInt16 uiSizeX, ezUInt16 uiSizeY, ezUInt16 uiSizeZ, const ezVec3* pSrc,
    ezVec3* pDst, ezVec3* pInputs, ezVec3* pScratch)
  {
    const ezUInt32 uiNumCells = (uiSizeX + 2) * (uiSizeY + 2) * (uiSizeZ + 2);

    CopyPreviousVelocity(pScratch, pSrc, pInputs, uiNumCells, fDampenFactor);

    Project3D(pScratch, pDst->GetData(), pDst->GetData() + uiNumCells, uiSizeX, uiSizeY, uiSizeZ);
    Advect3D(pDst, pScratch, uiSizeX, uiSizeY, uiSizeZ, deltaTime);
    Project3D(pDst, pScratch->GetData(), pScratch->GetData() + uiNumCells, uiSizeX, uiSizeY, uiSizeZ);
  }
} // namespace

//////////////////////////////////////////////////////////////////////////

ezWindSimulation::ezWindSimulation() = default;

ezWindSimulation::~ezWindSimulation()
{
  // make sure to finish the update task before destroying the object
  ezTaskSystem::WaitForGroup(m_UpdateTaskID);
}

void ezWindSimulation::Initialize(float fCellSize, ezUInt16 uiSizeX, ezUInt16 uiSizeY, ezUInt16 uiSizeZ /*= 1*/)
{
  m_fCellSize = fCellSize;
  m_fInverseCellSize = 1.0f / fCellSize;
  m_UpdateStep = ezTime::Milliseconds(100);

  m_uiSizeX = uiSizeX;
  m_uiSizeY = uiSizeY;
  m_uiSizeZ = uiSizeZ;

  m_uiIndexOffsetY = uiSizeX + 2;
  m_uiIndexOffsetZ = m_uiIndexOffsetY * (uiSizeY + 2);

  m_uiNumCells = 1;
  m_uiNumCells *= (m_uiSizeX + 2);
  m_uiNumCells *= (m_uiSizeY + 2);

  // Prev Velocity + Current Velocity + Next Velocity + Input Velocities + Scratch Data
  const ezUInt32 uiNumLayers = 5;

  if (IsVolumetric())
  {
    m_uiNumCells *= (m_uiSizeZ + 2);

    const ezUInt32 uiLayerSize = m_uiNumCells * 3;

    m_Values.SetCount(uiNumLayers * uiLayerSize);
    float* pCur = m_Values.GetData();

    m_pVelocities3D[0] = reinterpret_cast<ezVec3*>(pCur);
    pCur += uiLayerSize;
    m_pVelocities3D[1] = reinterpret_cast<ezVec3*>(pCur);
    pCur += uiLayerSize;
    m_pVelocities3D[2] = reinterpret_cast<ezVec3*>(pCur);
    pCur += uiLayerSize;
    m_pVelocityInputs3D = reinterpret_cast<ezVec3*>(pCur);
    pCur += uiLayerSize;
    m_pScratch3D = reinterpret_cast<ezVec3*>(pCur);
    pCur += uiLayerSize;
  }
  else
  {
    const ezUInt32 uiLayerSize = m_uiNumCells * 2;

    m_Values.SetCount(uiNumLayers * uiLayerSize);
    float* pCur = m_Values.GetData();

    m_pVelocities2D[0] = reinterpret_cast<ezVec2*>(pCur);
    pCur += uiLayerSize;
    m_pVelocities2D[1] = reinterpret_cast<ezVec2*>(pCur);
    pCur += uiLayerSize;
    m_pVelocities2D[2] = reinterpret_cast<ezVec2*>(pCur);
    pCur += uiLayerSize;
    m_pVelocityInputs2D = reinterpret_cast<ezVec2*>(pCur);
    pCur += uiLayerSize;
    m_pScratch2D = reinterpret_cast<ezVec2*>(pCur);
    pCur += uiLayerSize;
  }
}

void ezWindSimulation::Update(ezTime tDelta)
{
  EZ_PROFILE_SCOPE("Wind Simulation Update");

  const float fNewLerp = m_fUpdateFraction + tDelta.AsFloatInSeconds() / m_UpdateStep.AsFloatInSeconds();

  if (fNewLerp >= 1.0f)
  {
    // wait for the previous update task to finish (usually it should be done already)
    ezTaskSystem::WaitForGroup(m_UpdateTaskID);

    // update the current state indices
    m_fUpdateFraction = ezMath::Clamp(fNewLerp - 1.0f, 0.0f, 1.0f);
    m_uiPrevVelocities = m_uiCurVelocities;
    m_uiCurVelocities = m_uiNextVelocities;

    m_uiNextVelocities = (m_uiCurVelocities + 1) % 3;

    if (m_pUpdateTask == nullptr)
    {
      m_pUpdateTask = EZ_DEFAULT_NEW(ezDelegateTask<void>, "FluidWindUpdate", ezMakeDelegate(&ezWindSimulation::ComputeNextStep, this));
    }

    // launch a task to compute the next result in parallel
    // TODO: need more task priorities
    m_UpdateTaskID = ezTaskSystem::StartSingleTask(m_pUpdateTask.Borrow(), ezTaskPriority::LateNextFrame);
  }
  else
  {
    m_fUpdateFraction = fNewLerp;
  }
}

void ezWindSimulation::ComputeNextStep()
{
  const float deltaTime = m_UpdateStep.AsFloatInSeconds() / m_fCellSize;

  if (IsVolumetric())
  {
    StepWindSimulation3D(deltaTime, m_fDampenFactor, m_uiSizeX, m_uiSizeY, m_uiSizeZ, m_pVelocities3D[m_uiCurVelocities],
      m_pVelocities3D[m_uiNextVelocities], m_pVelocityInputs3D, m_pScratch3D);
  }
  else
  {
    StepWindSimulation2D(deltaTime, m_fDampenFactor, m_uiSizeX, m_uiSizeY, m_pVelocities2D[m_uiCurVelocities],
      m_pVelocities2D[m_uiNextVelocities], m_pVelocityInputs2D, m_pScratch2D);
  }
}

ezVec3 ezWindSimulation::MapPositionToCellIdx(const ezVec3& vPosition) const
{
  ezVec3 vCellIdx = vPosition * m_fInverseCellSize;
  vCellIdx.x += m_uiSizeX * 0.5f;
  vCellIdx.y += m_uiSizeY * 0.5f;
  vCellIdx.z += m_uiSizeZ * 0.5f;

  return vCellIdx;
}

ezVec2 ezWindSimulation::MapPositionToCellIdx(const ezVec2& vPosition) const
{
  ezVec2 vCellIdx = vPosition * m_fInverseCellSize;
  vCellIdx.x += m_uiSizeX * 0.5f;
  vCellIdx.y += m_uiSizeY * 0.5f;

  return vCellIdx;
}

ezVec2 ezWindSimulation::SampleVelocity2D(const ezVec2& vCellIdx) const
{
  const float minX = 0.5f;
  const float minY = 0.5f;
  const float maxX = m_uiSizeX + 0.5f;
  const float maxY = m_uiSizeY + 0.5f;

  // clamp sample position to valid range
  // allow border samples to enable fading result to zero
  const float x = ezMath::Clamp(vCellIdx.x, minX, maxX);
  const float y = ezMath::Clamp(vCellIdx.y, minY, maxY);

  // bilinear interpolate from the 4 sample position cells
  const ezUInt32 i0 = (ezUInt32)x;
  const ezUInt32 i1 = i0 + 1;
  const ezUInt32 j0 = (ezUInt32)y;
  const ezUInt32 j1 = j0 + 1;

  const float s1 = x - i0;
  const float s0 = 1 - s1;
  const float t1 = y - j0;
  const float t0 = 1 - t1;

  const ezUInt32 i0j0 = Idx(i0, j0);
  const ezUInt32 i0j1 = Idx(i0, j1);
  const ezUInt32 i1j0 = Idx(i1, j0);
  const ezUInt32 i1j1 = Idx(i1, j1);

  auto sampleComponent = [=](const ezVec2* pSrc) -> ezVec2 {
    return s0 * (t0 * pSrc[i0j0] + t1 * pSrc[i0j1]) + s1 * (t0 * pSrc[i1j0] + t1 * pSrc[i1j1]);
  };

  // probably not necessary at 10 Hz
  // return ezMath::Lerp(sampleComponent(m_pVelocities2D[m_uiPrevVelocities]), sampleComponent(m_pVelocities2D[m_uiCurVelocities]),
  // m_fLerpFactor);
  return sampleComponent(m_pVelocities2D[m_uiCurVelocities]);
}

ezVec3 ezWindSimulation::SampleVelocity3D(const ezVec3& vCellIdx) const
{
  const float minX = 0.5f;
  const float minY = 0.5f;
  const float minZ = 0.5f;
  const float maxX = m_uiSizeX + 0.5f;
  const float maxY = m_uiSizeY + 0.5f;
  const float maxZ = m_uiSizeZ + 0.5f;

  // clamp sample position to valid range
  // allow border samples to enable fading result to zero
  const float x = ezMath::Clamp(vCellIdx.x, minX, maxX);
  const float y = ezMath::Clamp(vCellIdx.y, minY, maxY);
  const float z = ezMath::Clamp(vCellIdx.z, minZ, maxZ);

  // trilinear interpolate from the 8 sample position cells
  const ezUInt32 i0 = (ezUInt32)x;
  const ezUInt32 i1 = i0 + 1;
  const ezUInt32 j0 = (ezUInt32)y;
  const ezUInt32 j1 = j0 + 1;
  const ezUInt32 k0 = (ezUInt32)z;
  const ezUInt32 k1 = k0 + 1;

  const float s1 = x - i0;
  const float s0 = 1 - s1;
  const float t1 = y - j0;
  const float t0 = 1 - t1;
  const float r1 = z - k0;
  const float r0 = 1 - r1;

  const ezUInt32 i0j0k0 = Idx(i0, j0, k0);
  const ezUInt32 i0j0k1 = Idx(i0, j0, k1);
  const ezUInt32 i0j1k0 = Idx(i0, j1, k0);
  const ezUInt32 i0j1k1 = Idx(i0, j1, k1);
  const ezUInt32 i1j0k0 = Idx(i1, j0, k0);
  const ezUInt32 i1j0k1 = Idx(i1, j0, k1);
  const ezUInt32 i1j1k0 = Idx(i1, j1, k0);
  const ezUInt32 i1j1k1 = Idx(i1, j1, k1);

  auto sampleComponent = [=](const ezVec3* pSrc) -> ezVec3 {
    return r0 * (s0 * (t0 * pSrc[i0j0k0] + t1 * pSrc[i0j1k0]) + s1 * (t0 * pSrc[i1j0k0] + t1 * pSrc[i1j1k0])) +
           r1 * (s0 * (t0 * pSrc[i0j0k1] + t1 * pSrc[i0j1k1]) + s1 * (t0 * pSrc[i1j0k1] + t1 * pSrc[i1j1k1]));
  };

  // probably not necessary at 10 Hz
  // return ezMath::Lerp(sampleComponent(m_pVelocities3D[m_uiPrevVelocities]), sampleComponent(m_pVelocities3D[m_uiCurVelocities]),
  // m_fLerpFactor);
  return sampleComponent(m_pVelocities3D[m_uiCurVelocities]);
}
