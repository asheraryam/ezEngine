#include <GameEnginePCH.h>

#include <GameEngine/Components/WindVolumeComponent.h>
#include <RendererCore/Debug/DebugRenderer.h>
#include <Foundation/Profiling/Profiling.h>

void ezWindVolumeComponent::addStream(int x, int y, float force)
{
  float* u = m_Wind.m_pVelocities[0];
  float* v = m_Wind.m_pVelocities[1];

  const int N = m_Wind.m_uiSizeX;

  u[m_Wind.Idx(x + 1, y)] = +force;
  v[m_Wind.Idx(x + 1, y)] = 0;

  u[m_Wind.Idx(x + 1, y - 1)] = +force * 0.7f;
  v[m_Wind.Idx(x + 1, y - 1)] = +force * 0.7f;
  u[m_Wind.Idx(x + 1, y + 1)] = +force * 0.7f;
  v[m_Wind.Idx(x + 1, y + 1)] = -force * 0.7f;
}

void ezWindVolumeComponent::addDrop(int x, int y, float force)
{
  float* u = m_Wind.m_pVelocities[0];
  float* v = m_Wind.m_pVelocities[1];

  const int N = m_Wind.m_uiSizeX;

  u[m_Wind.Idx(x - 1, y)] = -force;
  u[m_Wind.Idx(x + 1, y)] = +force;
  v[m_Wind.Idx(x, y - 1)] = -force;
  v[m_Wind.Idx(x, y + 1)] = +force;

  u[m_Wind.Idx(x - 1, y - 1)] = -force * 0.7f;
  v[m_Wind.Idx(x - 1, y - 1)] = +force * 0.7f;

  u[m_Wind.Idx(x + 1, y - 1)] = +force * 0.7f;
  v[m_Wind.Idx(x + 1, y - 1)] = +force * 0.7f;

  u[m_Wind.Idx(x + 1, y + 1)] = +force * 0.7f;
  v[m_Wind.Idx(x + 1, y + 1)] = -force * 0.7f;

  u[m_Wind.Idx(x - 1, y + 1)] = -force * 0.7f;
  v[m_Wind.Idx(x - 1, y + 1)] = -force * 0.7f;
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

  m_Wind.Initialize(128, 128);
}

void ezWindVolumeComponent::Update()
{
  m_Wind.Step();

  ezRandom rng;
  rng.InitializeFromCurrentTime();

  const int N = m_Wind.m_uiSizeX;

  {
    static int maxDrops = 20;
    if (rng.UIntInRange(100) < 1)
    {
      if (maxDrops > 0)
      {
        maxDrops--;
        addDrop(rng.UIntInRange(N-4) + 2, rng.UIntInRange(N-4) + 2, rng.FloatMinMax(1.5f, 10.0f));
      }
    }
  }

  {
    static int maxStream = 500;

    if (maxStream > 0)
    {
      maxStream--;
      addStream(10, 32, rng.FloatMinMax(0.05f, 0.1f));
    }
  }

  // draw_velocity
  {
    EZ_PROFILE_SCOPE("draw_velocity");

    int i, j;
    float x, y, h;

    h = 1.0f / N;

    ezDynamicArray<ezDebugRenderer::Line> lines;
    lines.Reserve(N * N);

    const float* u = m_Wind.m_pVelocities[0];
    const float* v = m_Wind.m_pVelocities[1];

    for (i = 1; i <= N; i++)
    {
      x = (i - 0.5f) * h;
      for (j = 1; j <= N; j++)
      {
        y = (j - 0.5f) * h;

        auto& l = lines.ExpandAndGetRef();
        l.m_start.Set(x, y, 0);
        l.m_end.Set(x + u[m_Wind.Idx(i, j)], y + v[m_Wind.Idx(i, j)], 0);
      }
    }

    ezDebugRenderer::DrawLines(GetWorld(), lines, ezColor::White);
  }
}
