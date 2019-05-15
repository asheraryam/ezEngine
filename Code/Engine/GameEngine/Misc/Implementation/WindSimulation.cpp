#include <GameEnginePCH.h>

#include <Foundation/Profiling/Profiling.h>
#include <GameEngine/Misc/WindSimulation.h>
#include <Foundation/Math/Vec3.h>

ezWindSimulation::ezWindSimulation() = default;
ezWindSimulation::~ezWindSimulation() = default;

#if 1

// uses the empty boundary to affect the simulation

#define ClampMinX 0
#define ClampMaxX (m_uiSizeX+1)
#define ClampMinY 0
#define ClampMaxY (m_uiSizeY+1)
#define ClampMinZ 0
#define ClampMaxZ (m_uiSizeZ+1)
#define ClampOffset 0.0f

#else

// clamps sampling points to the inner values (no boundary)

#define ClampMinX 1
#define ClampMaxX (m_uiSizeX)
#define ClampMinY 1
#define ClampMaxY (m_uiSizeY)
#define ClampMinZ 1
#define ClampMaxZ (m_uiSizeZ)
#define ClampOffset 1.0f

#endif


void ezWindSimulation::Initialize(float fCellSize, ezUInt16 uiSizeX, ezUInt16 uiSizeY, ezUInt16 uiSizeZ /*= 1*/)
{
  m_fCellSize = fCellSize;
  m_UpdateStep = ezTime::Milliseconds(100);

  m_uiSizeX = uiSizeX;
  m_uiSizeY = uiSizeY;
  m_uiSizeZ = uiSizeZ;

  m_uiIndexOffsetY = uiSizeX + 2;
  m_uiIndexOffsetZ = m_uiIndexOffsetY * (uiSizeY + 2);

  m_uiNumCells = 1;

  m_uiNumCells *= (m_uiSizeX + 2);
  m_uiNumCells *= (m_uiSizeY + 2);

  if (IsVolumetric())
  {
    m_uiNumCells *= (m_uiSizeZ + 2);
  }

  const ezUInt32 numFloats = m_uiNumCells * (IsVolumetric() ? 6 : 4);

  m_Values.SetCount(numFloats);

  float* pCur = m_Values.GetData();

  m_pVelocities[0] = pCur;
  pCur += m_uiNumCells;

  m_pVelocities[1] = pCur;
  pCur += m_uiNumCells;

  if (IsVolumetric())
  {
    m_pVelocities[2] = pCur;
    pCur += m_uiNumCells;
  }

  m_pPrevVelocities[0] = pCur;
  pCur += m_uiNumCells;

  m_pPrevVelocities[1] = pCur;
  pCur += m_uiNumCells;

  if (IsVolumetric())
  {
    m_pPrevVelocities[2] = pCur;
    pCur += m_uiNumCells;
  }
}

void ezWindSimulation::Step(ezTime tDelta)
{
  EZ_PROFILE_SCOPE("Wind Simulation");

  // TODO: use tDelta to advance internal interpolation factor

  CopyPreviousVelocity(m_pPrevVelocities[0], m_pVelocities[0]);
  CopyPreviousVelocity(m_pPrevVelocities[1], m_pVelocities[1]);

  if (IsVolumetric())
  {
    CopyPreviousVelocity(m_pPrevVelocities[2], m_pVelocities[2]);
    Project3D(m_pPrevVelocities[0], m_pPrevVelocities[1], m_pPrevVelocities[2], m_pVelocities[0], m_pVelocities[1]);
  }
  else
  {
    Project2D(m_pPrevVelocities[0], m_pPrevVelocities[1], m_pVelocities[0], m_pVelocities[1]);
  }

  Advect(m_pVelocities[0], m_pPrevVelocities[0]);
  Advect(m_pVelocities[1], m_pPrevVelocities[1]);

  if (IsVolumetric())
  {
    Advect(m_pVelocities[2], m_pPrevVelocities[2]);
    Project3D(m_pVelocities[0], m_pVelocities[1], m_pVelocities[2], m_pPrevVelocities[0], m_pPrevVelocities[1]);
  }
  else
  {
    Project2D(m_pVelocities[0], m_pVelocities[1], m_pPrevVelocities[0], m_pPrevVelocities[1]);
  }
}

void ezWindSimulation::CopyPreviousVelocity(float* pDst, const float* pSrc)
{
  const float fDampening = 0.995f;
  //const float fDampening = 1.0f;

  for (ezUInt32 i = 0; i < m_uiNumCells; ++i)
  {
    pDst[i] = pSrc[i] * fDampening;
  }
}

void ezWindSimulation::LinearSolve(float* pDst, const float* pPrev)
{
  const ezUInt32 uiNumIterations = 20;

  if (IsVolumetric())
  {
    for (ezUInt32 k = 0; k < uiNumIterations; k++)
    {
      for (ezInt32 z = 1; z <= m_uiSizeZ; ++z)
      {
        const ezInt32 zm = ezMath::Max<ezInt32>(z - 1, ClampMinZ);
        const ezInt32 zp = ezMath::Min<ezInt32>(z + 1, ClampMaxZ);

        for (ezInt32 y = 1; y <= m_uiSizeY; ++y)
        {
          const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
          const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

          for (ezInt32 x = 1; x <= m_uiSizeX; ++x)
          {
            const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
            const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

            pDst[Idx(x, y, z)] = (pPrev[Idx(x, y, z)] + (pDst[Idx(xm, y, z)] + pDst[Idx(xp, y, z)] + pDst[Idx(x, ym, z)] +
                                                          pDst[Idx(x, yp, z)] + pDst[Idx(x, y, zm)] + pDst[Idx(x, y, zp)])) /
                                 6.0f;
          }
        }
      }
    }
  }
  else
  {
    for (ezUInt32 k = 0; k < uiNumIterations; k++)
    {
      for (ezUInt32 y = 1; y <= m_uiSizeY; ++y)
      {
        const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
        const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

        for (ezUInt32 x = 1; x <= m_uiSizeX; ++x)
        {
          const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
          const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

          pDst[Idx(x, y, 0)] =
            (pPrev[Idx(x, y, 0)] + (pDst[Idx(xm, y, 0)] + pDst[Idx(xp, y, 0)] + pDst[Idx(x, ym, 0)] + pDst[Idx(x, yp, 0)])) * 0.25f;
        }
      }
    }
  }
}

void ezWindSimulation::Project2D(float* pDstU, float* pDstV, float* pScratch1, float* pScratch2)
{
  const float fNorm = -0.5f;

  for (ezUInt32 y = 1; y <= m_uiSizeY; ++y)
  {
    const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
    const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

    for (ezUInt32 x = 1; x <= m_uiSizeX; ++x)
    {
      const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
      const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

      pScratch2[Idx(x, y, 0)] = (pDstU[Idx(xp, y)] - pDstU[Idx(xm, y)] + pDstV[Idx(x, yp)] - pDstV[Idx(x, ym)]) * fNorm;
    }
  }

  ezMemoryUtils::ZeroFill(pScratch1, m_uiNumCells);

  LinearSolve(pScratch1, pScratch2);

  const float fNorm2 = 0.5f;

  for (ezUInt32 y = 1; y <= m_uiSizeY; ++y)
  {
    const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
    const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

    for (ezUInt32 x = 1; x <= m_uiSizeX; ++x)
    {
      const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
      const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

      pDstU[Idx(x, y)] -= fNorm2 * (pScratch1[Idx(xp, y)] - pScratch1[Idx(xm, y)]);
      pDstV[Idx(x, y)] -= fNorm2 * (pScratch1[Idx(x, yp)] - pScratch1[Idx(x, ym)]);
    }
  }
}

void ezWindSimulation::Project3D(float* pDstU, float* pDstV, float* pDstW, float* pScratch1, float* pScratch2)
{
  const float fNorm = -0.33f;

  for (ezUInt32 z = 1; z <= m_uiSizeZ; ++z)
  {
    const ezInt32 zm = ezMath::Max<ezInt32>(z - 1, ClampMinZ);
    const ezInt32 zp = ezMath::Min<ezInt32>(z + 1, ClampMaxZ);

    for (ezUInt32 y = 1; y <= m_uiSizeY; ++y)
    {
      const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
      const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

      for (ezUInt32 x = 1; x <= m_uiSizeX; ++x)
      {
        const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
        const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

        pScratch2[Idx(x, y, z)] = fNorm * (pDstU[Idx(xp, y, z)] - pDstU[Idx(xm, y, z)] + pDstV[Idx(x, yp, z)] - pDstV[Idx(x, ym, z)] +
          pDstW[Idx(x, y, zp)] - pDstW[Idx(x, y, zm)]);
      }
    }
  }

  ezMemoryUtils::ZeroFill(pScratch1, m_uiNumCells);

  LinearSolve(pScratch1, pScratch2);

  const float fNorm2 = 0.33f;

  for (ezUInt32 z = 1; z <= m_uiSizeZ; ++z)
  {
    const ezInt32 zm = ezMath::Max<ezInt32>(z - 1, ClampMinZ);
    const ezInt32 zp = ezMath::Min<ezInt32>(z + 1, ClampMaxZ);

    for (ezUInt32 y = 1; y <= m_uiSizeY; ++y)
    {
      const ezInt32 ym = ezMath::Max<ezInt32>(y - 1, ClampMinY);
      const ezInt32 yp = ezMath::Min<ezInt32>(y + 1, ClampMaxY);

      for (ezUInt32 x = 1; x <= m_uiSizeX; ++x)
      {
        const ezInt32 xm = ezMath::Max<ezInt32>(x - 1, ClampMinX);
        const ezInt32 xp = ezMath::Min<ezInt32>(x + 1, ClampMaxX);

        pDstU[Idx(x, y, z)] -= fNorm2 * (pScratch1[Idx(xp, y, z)] - pScratch1[Idx(xm, y, z)]);
        pDstV[Idx(x, y, z)] -= fNorm2 * (pScratch1[Idx(x, yp, z)] - pScratch1[Idx(x, ym, z)]);
        pDstW[Idx(x, y, z)] -= fNorm2 * (pScratch1[Idx(x, y, zp)] - pScratch1[Idx(x, y, zm)]);
      }
    }
  }
}

void ezWindSimulation::Advect(float* pDst, const float* pSrc)
{
  const float dt0 = m_UpdateStep.AsFloatInSeconds() / m_fCellSize;

  const float minX = 0.5f + ClampOffset;
  const float minY = 0.5f + ClampOffset;
  const float minZ = 0.5f + ClampOffset;
  const float maxX = m_uiSizeX + 0.5f - ClampOffset;
  const float maxY = m_uiSizeY + 0.5f - ClampOffset;
  const float maxZ = m_uiSizeZ + 0.5f - ClampOffset;

  if (IsVolumetric())
  {
    for (ezUInt32 k = 1; k <= m_uiSizeZ; ++k)
    {
      for (ezUInt32 j = 1; j <= m_uiSizeY; ++j)
      {
        for (ezUInt32 i = 1; i <= m_uiSizeX; ++i)
        {
          // compute reverse velocity sample position
          const float x = ezMath::Clamp(i - dt0 * m_pPrevVelocities[0][Idx(i, j, k)], minX, maxX);
          const float y = ezMath::Clamp(j - dt0 * m_pPrevVelocities[1][Idx(i, j, k)], minY, maxY);
          const float z = ezMath::Clamp(k - dt0 * m_pPrevVelocities[2][Idx(i, j, k)], minZ, maxZ);

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

          pDst[Idx(i, j, k)] = r0 * (s0 * (t0 * pSrc[Idx(i0, j0, k0)] + t1 * pSrc[Idx(i0, j1, k0)]) +
                                      s1 * (t0 * pSrc[Idx(i1, j0, k0)] + t1 * pSrc[Idx(i1, j1, k0)])) +
                               r1 * (s0 * (t0 * pSrc[Idx(i0, j0, k1)] + t1 * pSrc[Idx(i0, j1, k1)]) +
                                      s1 * (t0 * pSrc[Idx(i1, j0, k1)] + t1 * pSrc[Idx(i1, j1, k1)]));
        }
      }
    }
  }
  else
  {
    for (ezUInt32 j = 1; j <= m_uiSizeY; ++j)
    {
      for (ezUInt32 i = 1; i <= m_uiSizeX; ++i)
      {
        // compute reverse velocity sample position
        const float x = ezMath::Clamp(i - dt0 * m_pPrevVelocities[0][Idx(i, j)], minX, maxX);
        const float y = ezMath::Clamp(j - dt0 * m_pPrevVelocities[1][Idx(i, j)], minY, maxY);

        // bilinear interpolate from the 4 sample position cells
        const ezUInt32 i0 = (ezUInt32)x;
        const ezUInt32 i1 = i0 + 1;
        const ezUInt32 j0 = (ezUInt32)y;
        const ezUInt32 j1 = j0 + 1;

        const float s1 = x - i0;
        const float s0 = 1 - s1;
        const float t1 = y - j0;
        const float t0 = 1 - t1;

        pDst[Idx(i, j)] = s0 * (t0 * pSrc[Idx(i0, j0)] + t1 * pSrc[Idx(i0, j1)]) + s1 * (t0 * pSrc[Idx(i1, j0)] + t1 * pSrc[Idx(i1, j1)]);
      }
    }
  }
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

  auto sampleComponent = [=](const float* pSrc) -> float {
    return s0 * (t0 * pSrc[i0j0] + t1 * pSrc[i0j1]) + s1 * (t0 * pSrc[i1j0] + t1 * pSrc[i1j1]);
  };

  ezVec2 res;
  res.x = sampleComponent(m_pVelocities[0]);
  res.y = sampleComponent(m_pVelocities[1]);

  return res;
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

  auto sampleComponent = [=](const float* pSrc) -> float
  {
    return r0 * (s0 * (t0 * pSrc[i0j0k0] + t1 * pSrc[i0j1k0]) + s1 * (t0 * pSrc[i1j0k0] + t1 * pSrc[i1j1k0])) +
           r1 * (s0 * (t0 * pSrc[i0j0k1] + t1 * pSrc[i0j1k1]) + s1 * (t0 * pSrc[i1j0k1] + t1 * pSrc[i1j1k1]));
  };

  ezVec3 res;
  res.x = sampleComponent(m_pVelocities[0]);
  res.y = sampleComponent(m_pVelocities[1]);
  res.z = sampleComponent(m_pVelocities[2]);

  return res;  
}
