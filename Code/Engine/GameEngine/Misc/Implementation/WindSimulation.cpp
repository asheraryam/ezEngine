#include <GameEnginePCH.h>

#include <Foundation/Profiling/Profiling.h>
#include <GameEngine/Misc/WindSimulation.h>

ezWindSimulation::ezWindSimulation() = default;
ezWindSimulation::~ezWindSimulation() = default;

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

  AddTimeScaled(m_pPrevVelocities[0], m_pVelocities[0]);
  AddTimeScaled(m_pPrevVelocities[1], m_pVelocities[1]);

  if (IsVolumetric())
  {
    AddTimeScaled(m_pPrevVelocities[2], m_pVelocities[2]);
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

void ezWindSimulation::AddTimeScaled(float* pDst, const float* pSrc)
{
  const float fScale = m_UpdateStep.AsFloatInSeconds();
  const float fDampening = 1.0f; // 0.99f

  for (ezUInt32 i = 0; i < m_uiNumCells; ++i)
  {
    pDst[i] = fDampening * pSrc[i] + fScale * pDst[i];
  }

  ClearBounds(pDst);
}

void ezWindSimulation::ClearBounds(float* pDst)
{
  // 3D wall condition
  {
    // for (int j = 1; j < N - 1; j++)
    //{
    //  for (int i = 1; i < N - 1; i++)
    //  {
    //    x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
    //    x[IX(i, j, N - 1)] = b == 3 ? -x[IX(i, j, N - 2)] : x[IX(i, j, N - 2)];
    //  }
    //}
    // for (int k = 1; k < N - 1; k++)
    //{
    //  for (int i = 1; i < N - 1; i++)
    //  {
    //    x[IX(i, 0, k)] = b == 2 ? -x[IX(i, 1, k)] : x[IX(i, 1, k)];
    //    x[IX(i, N - 1, k)] = b == 2 ? -x[IX(i, N - 2, k)] : x[IX(i, N - 2, k)];
    //  }
    //}
    // for (int k = 1; k < N - 1; k++)
    //{
    //  for (int j = 1; j < N - 1; j++)
    //  {
    //    x[IX(0, j, k)] = b == 1 ? -x[IX(1, j, k)] : x[IX(1, j, k)];
    //    x[IX(N - 1, j, k)] = b == 1 ? -x[IX(N - 2, j, k)] : x[IX(N - 2, j, k)];
    //  }
    //}

    // x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, 1)]);
    // x[IX(0, N - 1, 0)] = 0.33f * (x[IX(1, N - 1, 0)] + x[IX(0, N - 2, 0)] + x[IX(0, N - 1, 1)]);
    // x[IX(0, 0, N - 1)] = 0.33f * (x[IX(1, 0, N - 1)] + x[IX(0, 1, N - 1)] + x[IX(0, 0, N)]);
    // x[IX(0, N - 1, N - 1)] = 0.33f * (x[IX(1, N - 1, N - 1)] + x[IX(0, N - 2, N - 1)] + x[IX(0, N - 1, N - 2)]);
    // x[IX(N - 1, 0, 0)] = 0.33f * (x[IX(N - 2, 0, 0)] + x[IX(N - 1, 1, 0)] + x[IX(N - 1, 0, 1)]);
    // x[IX(N - 1, N - 1, 0)] = 0.33f * (x[IX(N - 2, N - 1, 0)] + x[IX(N - 1, N - 2, 0)] + x[IX(N - 1, N - 1, 1)]);
    // x[IX(N - 1, 0, N - 1)] = 0.33f * (x[IX(N - 2, 0, N - 1)] + x[IX(N - 1, 1, N - 1)] + x[IX(N - 1, 0, N - 2)]);
    // x[IX(N - 1, N - 1, N - 1)] = 0.33f * (x[IX(N - 2, N - 1, N - 1)] + x[IX(N - 1, N - 2, N - 1)] + x[IX(N - 1, N - 1, N - 2)]);
  }

  // if (IsVolumetric())
  //{
  //  EZ_ASSERT_NOT_IMPLEMENTED;
  //}
  // else
  //{
  //  for (ezUInt16 i = 0; i < m_uiSizeX; ++i)
  //  {
  //    pDst[Idx(i, 0, 0)] = 0;
  //  }

  //  for (ezUInt16 i = 1; i < m_uiSizeY + 1; ++i)
  //  {
  //    pDst[Idx(0, i, 0)] = 0;
  //    pDst[Idx(m_uiSizeX + 1, i, 0)] = 0;
  //  }

  //  for (ezUInt16 i = 0; i < m_uiSizeX; ++i)
  //  {
  //    pDst[Idx(i, m_uiSizeY + 1, 0)] = 0;
  //  }
  //}
}

void ezWindSimulation::LinearSolve(float* pDst, const float* pPrev)
{
  const ezUInt32 uiNumIterations = 20;

  if (IsVolumetric())
  {
    for (ezUInt32 k = 0; k < uiNumIterations; k++)
    {
      for (ezUInt32 z = 1; z <= m_uiSizeZ; ++z)
      {
        for (ezUInt32 y = 1; y <= m_uiSizeY; ++y)
        {
          for (ezUInt32 x = 1; x <= m_uiSizeX; ++x)
          {
            pDst[Idx(x, y, z)] = (pPrev[Idx(x, y, z)] + (pDst[Idx(x - 1, y, z)] + pDst[Idx(x + 1, y, z)] + pDst[Idx(x, y - 1, z)] +
                                                          pDst[Idx(x, y + 1, z)] + pDst[Idx(x, y, z - 1)] + pDst[Idx(x, y, z + 1)])) /
                                 6.0f;
          }
        }
      }

      ClearBounds(pDst);
    }
  }
  else
  {
    for (ezUInt32 k = 0; k < uiNumIterations; k++)
    {
      for (ezUInt32 y = 1; y <= m_uiSizeY; ++y)
      {
        for (ezUInt32 x = 1; x <= m_uiSizeX; ++x)
        {
          pDst[Idx(x, y, 0)] =
            (pPrev[Idx(x, y, 0)] + (pDst[Idx(x - 1, y, 0)] + pDst[Idx(x + 1, y, 0)] + pDst[Idx(x, y - 1, 0)] + pDst[Idx(x, y + 1, 0)])) *
            0.25f;
        }
      }

      ClearBounds(pDst);
    }
  }
}

void ezWindSimulation::Project2D(float* pDstU, float* pDstV, float* pScratch1, float* pScratch2)
{
  const float fNorm = -0.5f;

  for (ezUInt32 y = 1; y <= m_uiSizeY; ++y)
  {
    for (ezUInt32 x = 1; x <= m_uiSizeX; ++x)
    {
      pScratch2[Idx(x, y, 0)] = (pDstU[Idx(x + 1, y)] - pDstU[Idx(x - 1, y)] + pDstV[Idx(x, y + 1)] - pDstV[Idx(x, y - 1)]) * fNorm;
    }
  }

  ClearBounds(pScratch2);

  ezMemoryUtils::ZeroFill(pScratch1, m_uiNumCells);

  LinearSolve(pScratch1, pScratch2);

  const float fNorm2 = 0.5f;

  for (ezUInt32 y = 1; y <= m_uiSizeY; ++y)
  {
    for (ezUInt32 x = 1; x <= m_uiSizeX; ++x)
    {
      pDstU[Idx(x, y)] -= fNorm2 * (pScratch1[Idx(x + 1, y)] - pScratch1[Idx(x - 1, y)]);
      pDstV[Idx(x, y)] -= fNorm2 * (pScratch1[Idx(x, y + 1)] - pScratch1[Idx(x, y - 1)]);
    }
  }

  ClearBounds(pDstU);
  ClearBounds(pDstV);
}

void ezWindSimulation::Project3D(float* pDstU, float* pDstV, float* pDstW, float* pScratch1, float* pScratch2)
{
  const float fNorm = -0.33f;

  for (ezUInt32 z = 1; z <= m_uiSizeZ; ++z)
  {
    for (ezUInt32 y = 1; y <= m_uiSizeY; ++y)
    {
      for (ezUInt32 x = 1; x <= m_uiSizeX; ++x)
      {
        pScratch2[Idx(x, y, z)] = fNorm * (pDstU[Idx(x + 1, y, z)] - pDstU[Idx(x - 1, y, z)] + pDstV[Idx(x, y + 1, z)] -
                                            pDstV[Idx(x, y - 1, z)] + pDstV[Idx(x, y, z + 1)] - pDstV[Idx(x, y, z - 1)]);
      }
    }
  }

  ClearBounds(pScratch2);

  ezMemoryUtils::ZeroFill(pScratch1, m_uiNumCells);

  LinearSolve(pScratch1, pScratch2);

  const float fNorm2 = 0.33f;

  for (ezUInt32 z = 1; z <= m_uiSizeZ; ++z)
  {
    for (ezUInt32 y = 1; y <= m_uiSizeY; ++y)
    {
      for (ezUInt32 x = 1; x <= m_uiSizeX; ++x)
      {
        pDstU[Idx(x, y, z)] -= fNorm2 * (pScratch1[Idx(x + 1, y, z)] - pScratch1[Idx(x - 1, y, z)]);
        pDstV[Idx(x, y, z)] -= fNorm2 * (pScratch1[Idx(x, y + 1, z)] - pScratch1[Idx(x, y - 1, z)]);
        pDstW[Idx(x, y, z)] -= fNorm2 * (pScratch1[Idx(x, y, z + 1)] - pScratch1[Idx(x, y, z - 1)]);
      }
    }
  }

  ClearBounds(pDstU);
  ClearBounds(pDstV);
  ClearBounds(pDstW);
}

void ezWindSimulation::Advect(float* pDst, const float* pSrc)
{
  const float dt0 = m_UpdateStep.AsFloatInSeconds() / m_fCellSize;

  const float maxX = m_uiSizeX + 0.5f;
  const float maxY = m_uiSizeY + 0.5f;
  const float maxZ = m_uiSizeZ + 0.5f;

  if (IsVolumetric())
  {
    for (ezUInt32 k = 1; k <= m_uiSizeZ; ++k)
    {
      for (ezUInt32 j = 1; j <= m_uiSizeY; ++j)
      {
        for (ezUInt32 i = 1; i <= m_uiSizeX; ++i)
        {
          // compute reverse velocity sample position
          const float x = ezMath::Clamp(i - dt0 * m_pPrevVelocities[0][Idx(i, j, k)], 0.5f, maxX);
          const float y = ezMath::Clamp(j - dt0 * m_pPrevVelocities[1][Idx(i, j, k)], 0.5f, maxY);
          const float z = ezMath::Clamp(k - dt0 * m_pPrevVelocities[2][Idx(i, j, k)], 0.5f, maxZ);

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
        const float x = ezMath::Clamp(i - dt0 * m_pPrevVelocities[0][Idx(i, j)], 0.5f, maxX);
        const float y = ezMath::Clamp(j - dt0 * m_pPrevVelocities[1][Idx(i, j)], 0.5f, maxY);

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

  ClearBounds(pDst);
}
