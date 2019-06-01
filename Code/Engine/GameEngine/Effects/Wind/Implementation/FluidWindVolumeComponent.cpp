#include <GameEnginePCH.h>

#include <Foundation/Profiling/Profiling.h>
#include <GameEngine/Effects/Wind/FluidWindVolumeComponent.h>
#include <RendererCore/Debug/DebugRenderer.h>

#include <Core/WorldSerializer/WorldReader.h>
#include <Core/WorldSerializer/WorldWriter.h>

// clang-format off
EZ_BEGIN_COMPONENT_TYPE(ezFluidWindVolumeComponent, 1, ezComponentMode::Static)
{
  EZ_BEGIN_PROPERTIES
  {
    EZ_MEMBER_PROPERTY("CellsX", m_uiCellsX)->AddAttributes(new ezDefaultValueAttribute(32), new ezClampValueAttribute(1, 255)),
    EZ_MEMBER_PROPERTY("CellsY", m_uiCellsY)->AddAttributes(new ezDefaultValueAttribute(32), new ezClampValueAttribute(1, 255)),
    EZ_MEMBER_PROPERTY("CellsZ", m_uiCellsZ)->AddAttributes(new ezDefaultValueAttribute(4), new ezClampValueAttribute(1, 255)),
    EZ_MEMBER_PROPERTY("CellsSize", m_fCellSize)->AddAttributes(new ezDefaultValueAttribute(0.5f)),
    EZ_MEMBER_PROPERTY("Visualize", m_bVisualize),
  }
  EZ_END_PROPERTIES;
  EZ_BEGIN_ATTRIBUTES
  {
    new ezCategoryAttribute("Effects/Wind"),
  }
  EZ_END_ATTRIBUTES;
}
EZ_END_DYNAMIC_REFLECTED_TYPE;
// clang-format on

ezFluidWindVolumeComponent::ezFluidWindVolumeComponent() = default;
ezFluidWindVolumeComponent::~ezFluidWindVolumeComponent() = default;

void ezFluidWindVolumeComponent::SerializeComponent(ezWorldWriter& stream) const
{
  SUPER::SerializeComponent(stream);
  auto& s = stream.GetStream();

  // Version 1
  s << m_uiCellsX;
  s << m_uiCellsY;
  s << m_uiCellsZ;
  s << m_fCellSize;
  s << m_bVisualize;
}

void ezFluidWindVolumeComponent::DeserializeComponent(ezWorldReader& stream)
{
  SUPER::DeserializeComponent(stream);
  // const ezUInt32 uiVersion = stream.GetComponentTypeVersion(GetStaticRTTI());
  auto& s = stream.GetStream();

  // Version 1
  s >> m_uiCellsX;
  s >> m_uiCellsY;
  s >> m_uiCellsZ;
  s >> m_fCellSize;
  s >> m_bVisualize;
}

static int maxDrops = 20;
static int maxStream = 500;

void ezFluidWindVolumeComponent::OnSimulationStarted()
{
  EZ_ASSERT_DEV(m_pWindModule == nullptr, "Wind module should be null");

  maxDrops = 20;
  maxStream = 500;

  m_FluidVolume.m_Simulation.Initialize(m_fCellSize, m_uiCellsX, m_uiCellsY, m_uiCellsZ);

  m_pWindModule = GetWorld()->GetOrCreateModule<ezWindWorldModuleInterface>();
  m_pWindModule->AddWindVolume(&m_FluidVolume);
}

void ezFluidWindVolumeComponent::OnDeactivated()
{
  if (m_pWindModule != nullptr)
  {
    m_pWindModule->RemoveWindVolume(&m_FluidVolume);
    m_pWindModule = nullptr;
  }
}

void ezFluidWindVolumeComponent::Update()
{
  const ezTransform ownTransform = GetOwner()->GetGlobalTransform();

  m_FluidVolume.m_Simulation.Update(GetWorld()->GetClock().GetTimeDiff());
  m_FluidVolume.m_vPosition = ownTransform.m_vPosition;
  m_FluidVolume.m_qRotation = ownTransform.m_qRotation;

  if (m_bVisualize)
  {
    EZ_PROFILE_SCOPE("draw_velocity");

    ezDynamicArray<ezDebugRenderer::Line> lines;
    lines.Reserve(m_FluidVolume.m_Simulation.GetSizeX() * m_FluidVolume.m_Simulation.GetSizeY() * m_FluidVolume.m_Simulation.GetSizeZ());


    const float h = m_FluidVolume.m_Simulation.GetCellSize();
    const float h2 = h * 0.5f;

    ezVec3 vOffset;
    vOffset.x = (m_FluidVolume.m_Simulation.GetSizeX() * -0.5f) * m_FluidVolume.m_Simulation.GetCellSize();
    vOffset.y = (m_FluidVolume.m_Simulation.GetSizeY() * -0.5f) * m_FluidVolume.m_Simulation.GetCellSize();
    vOffset.z = (m_FluidVolume.m_Simulation.GetSizeZ() * -0.5f) * m_FluidVolume.m_Simulation.GetCellSize();

    const ezTransform ownerTransform = GetOwner()->GetGlobalTransform();

    if (m_FluidVolume.m_Simulation.IsVolumetric())
    {
      const ezVec3* vel = m_FluidVolume.m_Simulation.GetVelocities3D();

      for (int k = 1; k <= m_FluidVolume.m_Simulation.GetSizeZ(); k++)
      {
        const float z = (k - 0.5f) * h;

        for (int j = 1; j <= m_FluidVolume.m_Simulation.GetSizeY(); j++)
        {
          const float y = (j - 0.5f) * h;

          for (int i = 1; i <= m_FluidVolume.m_Simulation.GetSizeX(); i++)
          {
            const float x = (i - 0.5f) * h;

            auto& l = lines.ExpandAndGetRef();
            l.m_start.Set(x, y, z);
            l.m_start += vOffset;
            l.m_start = ownerTransform * l.m_start;

            const ezVec3 sample = vel[m_FluidVolume.m_Simulation.Idx(i, j, k)];

            l.m_end.Set(x + sample.x * h2, y + sample.y * h2, z + sample.z * h2);

            l.m_end += vOffset;
            l.m_end = ownerTransform * l.m_end;
          }
        }
      }
    }
    else
    {
      const ezVec2* vel = m_FluidVolume.m_Simulation.GetVelocities2D();

      for (int j = 1; j <= m_FluidVolume.m_Simulation.GetSizeY(); j++)
      {
        const float y = (j - 0.5f) * h;

        for (int i = 1; i <= m_FluidVolume.m_Simulation.GetSizeX(); i++)
        {
          const float x = (i - 0.5f) * h;

          auto& l = lines.ExpandAndGetRef();
          l.m_start.Set(x, y, 0);
          l.m_start += vOffset;
          l.m_start = ownerTransform * l.m_start;

          const ezVec2 sample = vel[m_FluidVolume.m_Simulation.Idx(i, j)];

          l.m_end.Set(x + sample.x * h2, y + sample.y * h2, 0);

          l.m_end += vOffset;
          l.m_end = ownerTransform * l.m_end;
        }
      }
    }

    ezDebugRenderer::DrawLines(GetWorld(), lines, ezColor::White);
  }
}

//////////////////////////////////////////////////////////////////////////

ezResult ezFluidWindVolume::GetWindAt(const ezVec3& vGlobalPosition, ezVec3& out_vWind)
{
  const ezVec3 vLocalPos = -m_qRotation * (vGlobalPosition - m_vPosition);
  const ezVec3 vCellIdx = m_Simulation.MapPositionToCellIdx(vLocalPos);

  if (m_Simulation.IsVolumetric())
    out_vWind = m_Simulation.SampleVelocity3D(vCellIdx);
  else
    out_vWind = m_Simulation.SampleVelocity2D(vCellIdx.GetAsVec2()).GetAsVec3(0.0f);

  out_vWind *= 10.0f;

  return EZ_SUCCESS;
}

void ezFluidWindVolume::ApplyForceSphere(const ezVec3& vCenter0, float fRadius, float fStrength)
{
  const ezVec3 vCenter = vCenter0 - m_vPosition;

  const float fCellSize = m_Simulation.GetCellSize();
  const ezInt32 iRadius = ezMath::Ceil(fRadius / fCellSize);
  const ezInt32 iSize = iRadius * 2;

  ezDynamicArray<ezVec3> samplePoints;
  samplePoints.Reserve(iSize * iSize * iSize * 3 / 4);

  const float fRadiusSqr = ezMath::Square(fRadius);

  for (ezInt32 z = -iRadius; z <= +iRadius; ++z)
  {
    for (ezInt32 y = -iRadius; y <= +iRadius; ++y)
    {
      for (ezInt32 x = -iRadius; x <= +iRadius; ++x)
      {
        const ezVec3 pos(x * fCellSize, y * fCellSize, z * fCellSize);

        if (pos.GetLengthSquared() < fRadiusSqr)
          samplePoints.PushBack(pos);
      }
    }
  }

  const ezVec3 vLocalCenter = m_Simulation.MapPositionToCellIdx(vCenter);


  for (const ezVec3 pos : samplePoints)
  {
    ezVec3 finalPos = vCenter + pos;

    const ezVec3 cellIdx = m_Simulation.MapPositionToCellIdx(finalPos);
    ezVec3 snappedPos;
    snappedPos.x = ezMath::Floor(cellIdx.x);
    snappedPos.y = ezMath::Floor(cellIdx.y);
    snappedPos.z = ezMath::Floor(cellIdx.z);

    if (snappedPos.x < 0 || snappedPos.x >= m_Simulation.GetSizeX())
      continue;
    if (snappedPos.y < 0 || snappedPos.y >= m_Simulation.GetSizeY())
      continue;
    if (snappedPos.z < 0 || snappedPos.z >= m_Simulation.GetSizeZ())
      continue;

    ezVec3 dir = snappedPos - vLocalCenter;
    if (dir.NormalizeIfNotZero(ezVec3::ZeroVector()).Failed())
      continue;

    dir *= fStrength;

    ezUInt32 idx = m_Simulation.Idx(snappedPos.x, snappedPos.y, snappedPos.z);
    ezVec3& vel = m_Simulation.GetVelocityInputs3D()[idx];
    vel = dir;
  }
}
