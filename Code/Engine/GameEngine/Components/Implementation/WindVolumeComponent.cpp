#include <GameEnginePCH.h>

#include <Foundation/Profiling/Profiling.h>
#include <GameEngine/Components/WindVolumeComponent.h>
#include <RendererCore/Debug/DebugRenderer.h>

void ezWindVolumeComponent::addStream(int x, int y, float force)
{
  float* u = m_Wind.GetVelocitiesX();
  float* v = m_Wind.GetVelocitiesY();
  const int z = m_Wind.IsVolumetric() ? m_Wind.GetSizeZ() / 2 : 0;

  u[m_Wind.Idx(x + 1, y, z)] = +force;
  v[m_Wind.Idx(x + 1, y, z)] = 0;

  u[m_Wind.Idx(x + 1, y - 1, z)] = +force * 0.7f;
  v[m_Wind.Idx(x + 1, y - 1, z)] = +force * 0.7f;
  u[m_Wind.Idx(x + 1, y + 1, z)] = +force * 0.7f;
  v[m_Wind.Idx(x + 1, y + 1, z)] = -force * 0.7f;
}

void ezWindVolumeComponent::addDrop(int x, int y, float force)
{
  float* u = m_Wind.GetVelocitiesX();
  float* v = m_Wind.GetVelocitiesY();
  const int z = m_Wind.IsVolumetric() ? m_Wind.GetSizeZ() / 2 : 0;

  u[m_Wind.Idx(x - 1, y, z)] = -force;
  u[m_Wind.Idx(x + 1, y, z)] = +force;
  v[m_Wind.Idx(x, y - 1, z)] = -force;
  v[m_Wind.Idx(x, y + 1, z)] = +force;

  u[m_Wind.Idx(x - 1, y - 1, z)] = -force * 0.7f;
  v[m_Wind.Idx(x - 1, y - 1, z)] = +force * 0.7f;

  u[m_Wind.Idx(x + 1, y - 1, z)] = +force * 0.7f;
  v[m_Wind.Idx(x + 1, y - 1, z)] = +force * 0.7f;

  u[m_Wind.Idx(x + 1, y + 1, z)] = +force * 0.7f;
  v[m_Wind.Idx(x + 1, y + 1, z)] = -force * 0.7f;

  u[m_Wind.Idx(x - 1, y + 1, z)] = -force * 0.7f;
  v[m_Wind.Idx(x - 1, y + 1, z)] = -force * 0.7f;
}


//////////////////////////////////////////////////////////////////////////

#include <Core/WorldSerializer/WorldReader.h>
#include <Core/WorldSerializer/WorldWriter.h>

// clang-format off
EZ_BEGIN_COMPONENT_TYPE(ezWindVolumeComponent, 1, ezComponentMode::Dynamic)
{
  //EZ_BEGIN_PROPERTIES
  //{
  //  EZ_ACCESSOR_PROPERTY("Gradient", GetColorGradientFile, SetColorGradientFile)->AddAttributes(new ezAssetBrowserAttribute("ColorGradient")),
  //  EZ_MEMBER_PROPERTY("Duration", m_Duration),
  //}
  //EZ_END_PROPERTIES;
  //EZ_BEGIN_ATTRIBUTES
  //{
  //  new ezCategoryAttribute("Animation"),
  //}
  //EZ_END_ATTRIBUTES;
}
EZ_END_DYNAMIC_REFLECTED_TYPE;
// clang-format on

ezWindVolumeComponent::ezWindVolumeComponent() {}

void ezWindVolumeComponent::SerializeComponent(ezWorldWriter& stream) const
{
  SUPER::SerializeComponent(stream);
  auto& s = stream.GetStream();
}

void ezWindVolumeComponent::DeserializeComponent(ezWorldReader& stream)
{
  SUPER::DeserializeComponent(stream);
  // const ezUInt32 uiVersion = stream.GetComponentTypeVersion(GetStaticRTTI());
  auto& s = stream.GetStream();
}

void ezWindVolumeComponent::OnSimulationStarted()
{
  SUPER::OnSimulationStarted();

  m_Wind.Initialize(0.25f, 64, 64, 8);
}

void ezWindVolumeComponent::Update()
{
  ezRandom rng;
  rng.InitializeFromCurrentTime();

  {
    static int maxDrops = 20;
    if (rng.UIntInRange(100) < 1)
    {
      if (maxDrops > 0)
      {
        maxDrops--;
        //addDrop(rng.UIntInRange(m_Wind.GetSizeX() - 4) + 2, rng.UIntInRange(m_Wind.GetSizeY() - 4) + 2, rng.FloatMinMax(m_Wind.GetCellSize() * 50.0f, m_Wind.GetCellSize() * 100.0f));
      }
    }
  }

  {
    static int maxStream = 500;

    if (maxStream > 0)
    {
      maxStream--;
      addStream(10, 32, rng.FloatMinMax(m_Wind.GetCellSize() * 1000.0f, m_Wind.GetCellSize() * 1100.0f));
    }
  }

  m_Wind.Step(ezTime::Milliseconds(100));

  // draw_velocity
  {
    EZ_PROFILE_SCOPE("draw_velocity");

    ezDynamicArray<ezDebugRenderer::Line> lines;
    lines.Reserve(m_Wind.GetSizeX() * m_Wind.GetSizeY() * m_Wind.GetSizeZ());

    const float* u = m_Wind.GetVelocitiesX();
    const float* v = m_Wind.GetVelocitiesY();
    const float* w = m_Wind.GetVelocitiesZ();

    const float h = m_Wind.GetCellSize();
    const float h2 = h * 0.5f;

    for (int k0 = 1; k0 <= m_Wind.GetSizeZ(); k0++)
    {
      const int k = m_Wind.IsVolumetric() ? k0 : 0;
      const float z = (k - 0.5f) * h;

      for (int j = 1; j <= m_Wind.GetSizeY(); j++)
      {
        const float y = (j - 0.5f) * h;

        for (int i = 1; i <= m_Wind.GetSizeX(); i++)
        {
          const float x = (i - 0.5f) * h;

          auto& l = lines.ExpandAndGetRef();
          l.m_start.Set(x, y, z);

          if (w != nullptr)
          {
            l.m_end.Set(x + u[m_Wind.Idx(i, j, k)] * h2, y + v[m_Wind.Idx(i, j, k)] * h2, z + w[m_Wind.Idx(i, j, k)] * h2);
          }
          else
          {
            l.m_end.Set(x + u[m_Wind.Idx(i, j)] * h2, y + v[m_Wind.Idx(i, j)] * h2, z);
          }
        }
      }
    }

    ezDebugRenderer::DrawLines(GetWorld(), lines, ezColor::White);
  }
}
