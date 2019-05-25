#include <GameEnginePCH.h>

#include <Foundation/Profiling/Profiling.h>
#include <GameEngine/Components/WindVolumeComponent.h>
#include <RendererCore/Debug/DebugRenderer.h>

void ezWindVolumeComponent::addStream(int x, int y, float force)
{
  float* u = m_FluidVolume.m_Simulation.GetVelocitiesX();
  float* v = m_FluidVolume.m_Simulation.GetVelocitiesY();
  const int z = m_FluidVolume.m_Simulation.IsVolumetric() ? m_FluidVolume.m_Simulation.GetSizeZ() / 2 : 0;

  u[m_FluidVolume.m_Simulation.Idx(x + 1, y, z)] = +force;
  v[m_FluidVolume.m_Simulation.Idx(x + 1, y, z)] = 0;

  u[m_FluidVolume.m_Simulation.Idx(x + 1, y - 1, z)] = +force * 0.7f;
  v[m_FluidVolume.m_Simulation.Idx(x + 1, y - 1, z)] = +force * 0.7f;
  u[m_FluidVolume.m_Simulation.Idx(x + 1, y + 1, z)] = +force * 0.7f;
  v[m_FluidVolume.m_Simulation.Idx(x + 1, y + 1, z)] = -force * 0.7f;
}

void ezWindVolumeComponent::addDrop(int x, int y, float force)
{
  float* u = m_FluidVolume.m_Simulation.GetVelocitiesX();
  float* v = m_FluidVolume.m_Simulation.GetVelocitiesY();
  const int z = m_FluidVolume.m_Simulation.IsVolumetric() ? m_FluidVolume.m_Simulation.GetSizeZ() / 2 : 0;

  u[m_FluidVolume.m_Simulation.Idx(x - 1, y, z)] = -force;
  u[m_FluidVolume.m_Simulation.Idx(x + 1, y, z)] = +force;
  v[m_FluidVolume.m_Simulation.Idx(x, y - 1, z)] = -force;
  v[m_FluidVolume.m_Simulation.Idx(x, y + 1, z)] = +force;

  u[m_FluidVolume.m_Simulation.Idx(x - 1, y - 1, z)] = -force * 0.7f;
  v[m_FluidVolume.m_Simulation.Idx(x - 1, y - 1, z)] = +force * 0.7f;

  u[m_FluidVolume.m_Simulation.Idx(x + 1, y - 1, z)] = +force * 0.7f;
  v[m_FluidVolume.m_Simulation.Idx(x + 1, y - 1, z)] = +force * 0.7f;

  u[m_FluidVolume.m_Simulation.Idx(x + 1, y + 1, z)] = +force * 0.7f;
  v[m_FluidVolume.m_Simulation.Idx(x + 1, y + 1, z)] = -force * 0.7f;

  u[m_FluidVolume.m_Simulation.Idx(x - 1, y + 1, z)] = -force * 0.7f;
  v[m_FluidVolume.m_Simulation.Idx(x - 1, y + 1, z)] = -force * 0.7f;
}


//////////////////////////////////////////////////////////////////////////

#include <Core/WorldSerializer/WorldReader.h>
#include <Core/WorldSerializer/WorldWriter.h>

// clang-format off
EZ_BEGIN_COMPONENT_TYPE(ezWindVolumeComponent, 1, ezComponentMode::Static)
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
    new ezCategoryAttribute("Effects"),
  }
  EZ_END_ATTRIBUTES;
}
EZ_END_DYNAMIC_REFLECTED_TYPE;
// clang-format on

ezWindVolumeComponent::ezWindVolumeComponent() = default;
ezWindVolumeComponent::~ezWindVolumeComponent() = default;

void ezWindVolumeComponent::SerializeComponent(ezWorldWriter& stream) const
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

void ezWindVolumeComponent::DeserializeComponent(ezWorldReader& stream)
{
  SUPER::DeserializeComponent(stream);
  //const ezUInt32 uiVersion = stream.GetComponentTypeVersion(GetStaticRTTI());
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

void ezWindVolumeComponent::OnActivated()
{
  SUPER::OnSimulationStarted();

  maxDrops = 20;
  maxStream = 500;

  m_FluidVolume.m_Simulation.Initialize(m_fCellSize, m_uiCellsX, m_uiCellsY, m_uiCellsZ);

  m_pWindModule = GetWorld()->GetOrCreateModule<ezWindWorldModuleInterface>();
  m_pWindModule->AddWindVolume(&m_FluidVolume);
}

void ezWindVolumeComponent::OnDeactivated()
{
  if (m_pWindModule != nullptr)
  {
    m_pWindModule->RemoveWindVolume(&m_FluidVolume);
    m_pWindModule = nullptr;
  }
}

void ezWindVolumeComponent::Update()
{
  ezRandom rng;
  rng.InitializeFromCurrentTime();

  {
    if (rng.UIntInRange(100) < 1)
    {
      if (maxDrops > 0)
      {
        maxDrops--;
        addDrop(rng.UIntInRange(m_FluidVolume.m_Simulation.GetSizeX() - 4) + 2,
          rng.UIntInRange(m_FluidVolume.m_Simulation.GetSizeY() - 4) + 2,
          rng.FloatMinMax(m_FluidVolume.m_Simulation.GetCellSize() * 1000.0f, m_FluidVolume.m_Simulation.GetCellSize() * 10000.0f));
      }
    }
  }

  {

    if (maxStream > 0)
    {
      maxStream--;
      addStream(
        10, 32, rng.FloatMinMax(m_FluidVolume.m_Simulation.GetCellSize() * 20.0f, m_FluidVolume.m_Simulation.GetCellSize() * 100.0f));
    }
  }

  const ezTransform ownTransform = GetOwner()->GetGlobalTransform();

  m_FluidVolume.m_Simulation.Step(ezTime::Milliseconds(100));
  m_FluidVolume.m_vPosition = ownTransform.m_vPosition;
  m_FluidVolume.m_qRotation = ownTransform.m_qRotation;

  // draw_velocity
  if (m_bVisualize)
  {
    EZ_PROFILE_SCOPE("draw_velocity");

    ezDynamicArray<ezDebugRenderer::Line> lines;
    lines.Reserve(m_FluidVolume.m_Simulation.GetSizeX() * m_FluidVolume.m_Simulation.GetSizeY() * m_FluidVolume.m_Simulation.GetSizeZ());

    const float* u = m_FluidVolume.m_Simulation.GetVelocitiesX();
    const float* v = m_FluidVolume.m_Simulation.GetVelocitiesY();
    const float* w = m_FluidVolume.m_Simulation.GetVelocitiesZ();

    const float h = m_FluidVolume.m_Simulation.GetCellSize();
    const float h2 = h * 0.5f;

    ezVec3 vOffset;
    vOffset.x = (m_FluidVolume.m_Simulation.GetSizeX() * -0.5f) * m_FluidVolume.m_Simulation.GetCellSize();
    vOffset.y = (m_FluidVolume.m_Simulation.GetSizeY() * -0.5f) * m_FluidVolume.m_Simulation.GetCellSize();
    vOffset.z = (m_FluidVolume.m_Simulation.GetSizeZ() * -0.5f) * m_FluidVolume.m_Simulation.GetCellSize();

    const ezTransform ownerTransform = GetOwner()->GetGlobalTransform();

    for (int k0 = 1; k0 <= m_FluidVolume.m_Simulation.GetSizeZ(); k0++)
    {
      const int k = m_FluidVolume.m_Simulation.IsVolumetric() ? k0 : 0;
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

          if (w != nullptr)
          {
            l.m_end.Set(x + u[m_FluidVolume.m_Simulation.Idx(i, j, k)] * h2, y + v[m_FluidVolume.m_Simulation.Idx(i, j, k)] * h2,
              z + w[m_FluidVolume.m_Simulation.Idx(i, j, k)] * h2);
          }
          else
          {
            l.m_end.Set(x + u[m_FluidVolume.m_Simulation.Idx(i, j)] * h2, y + v[m_FluidVolume.m_Simulation.Idx(i, j)] * h2, z);
          }

          l.m_end += vOffset;
          l.m_end = ownerTransform * l.m_end;
        }
      }
    }

    ezDebugRenderer::DrawLines(GetWorld(), lines, ezColor::White);
  }
}

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
