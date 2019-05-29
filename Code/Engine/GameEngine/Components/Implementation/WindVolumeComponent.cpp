#include <GameEnginePCH.h>

#include <Foundation/Profiling/Profiling.h>
#include <GameEngine/Components/WindVolumeComponent.h>
#include <RendererCore/Debug/DebugRenderer.h>

#define AddForce(cell, force) cell = ezMath::Max(cell, force)

void ezWindVolumeComponent::addStream(int x, int y, float force)
{
  if (m_FluidVolume.m_Simulation.IsVolumetric())
  {
    ezVec3* vel = m_FluidVolume.m_Simulation.GetVelocityInputs3D();
    const int z = m_FluidVolume.m_Simulation.GetSizeZ() / 2;

    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y, z)], ezVec3(force, 0, 0));

    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y - 1, z)], ezVec3(force * 0.7f, force * 0.7f, 0));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y + 1, z)], ezVec3(force * 0.7f, -force * 0.7f, 0));
  }
  else
  {
    ezVec2* vel = m_FluidVolume.m_Simulation.GetVelocityInputs2D();

    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y)], ezVec2(force, 0));

    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y - 1)], ezVec2(force * 0.7f, force * 0.7f));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y + 1)], ezVec2(force * 0.7f, -force * 0.7f));
  }
}

void ezWindVolumeComponent::addDrop(int x, int y, float force)
{
  if (m_FluidVolume.m_Simulation.IsVolumetric())
  {
    ezVec3* vel = m_FluidVolume.m_Simulation.GetVelocityInputs3D();
    const int z = m_FluidVolume.m_Simulation.GetSizeZ() / 2;

    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x - 1, y, z)], ezVec3(-force, 0, 0));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y, z)], ezVec3(+force, 0, 0));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x, y - 1, z)], ezVec3(0, -force, 0));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x, y + 1, z)], ezVec3(0, +force, 0));

    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x - 1, y - 1, z)], ezVec3(-force * 0.7f, +force * 0.7f, 0));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y - 1, z)], ezVec3(+force * 0.7f, +force * 0.7f, 0));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y + 1, z)], ezVec3(+force * 0.7f, -force * 0.7f, 0));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x - 1, y + 1, z)], ezVec3(-force * 0.7f, -force * 0.7f, 0));
  }
  else
  {
    ezVec2* vel = m_FluidVolume.m_Simulation.GetVelocityInputs2D();

    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x - 1, y)], ezVec2(-force, 0));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y)], ezVec2(+force, 0));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x, y - 1)], ezVec2(0, -force));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x, y + 1)], ezVec2(0, +force));

    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x - 1, y - 1)], ezVec2(-force * 0.7f, +force * 0.7f));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y - 1)], ezVec2(+force * 0.7f, +force * 0.7f));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x + 1, y + 1)], ezVec2(+force * 0.7f, -force * 0.7f));
    AddForce(vel[m_FluidVolume.m_Simulation.Idx(x - 1, y + 1)], ezVec2(-force * 0.7f, -force * 0.7f));
  }
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

void ezWindVolumeComponent::OnSimulationStarted()
{
  EZ_ASSERT_DEV(m_pWindModule == nullptr, "Wind module should be null");

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
