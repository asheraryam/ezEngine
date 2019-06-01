#include <GameEnginePCH.h>

#include <Core/WorldSerializer/WorldReader.h>
#include <Core/WorldSerializer/WorldWriter.h>
#include <GameEngine/Effects/Wind/WindForceComponent.h>
#include <RendererCore/Debug/DebugRenderer.h>

// clang-format off
EZ_BEGIN_COMPONENT_TYPE(ezWindForceComponent, 1, ezComponentMode::Static)
{
  EZ_BEGIN_PROPERTIES
  {
    EZ_MEMBER_PROPERTY("Radius", m_fRadius)->AddAttributes(new ezDefaultValueAttribute(1.0f), new ezClampValueAttribute(0.0f, 1000.0f)),
    EZ_MEMBER_PROPERTY("Strength", m_fStrength)->AddAttributes(new ezDefaultValueAttribute(1.0f), new ezClampValueAttribute(0.0f, 10000.0f)),
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

ezWindForceComponent::ezWindForceComponent() = default;
ezWindForceComponent::~ezWindForceComponent() = default;

void ezWindForceComponent::SerializeComponent(ezWorldWriter& stream) const
{
  SUPER::SerializeComponent(stream);
  auto& s = stream.GetStream();

  // Version 1
  s << m_fRadius;
  s << m_fStrength;
}

void ezWindForceComponent::DeserializeComponent(ezWorldReader& stream)
{
  SUPER::DeserializeComponent(stream);
  // const ezUInt32 uiVersion = stream.GetComponentTypeVersion(GetStaticRTTI());
  auto& s = stream.GetStream();

  // Version 1
  s >> m_fRadius;
  s >> m_fStrength;
}

void ezWindForceComponent::OnSimulationStarted()
{
  EZ_ASSERT_DEV(m_pWindModule == nullptr, "Wind module should be null");
  m_pWindModule = GetWorld()->GetOrCreateModule<ezWindWorldModuleInterface>();
}

void ezWindForceComponent::Update()
{
  m_pWindModule->ApplyForceSphere(GetOwner()->GetGlobalPosition(), m_fRadius, m_fStrength);
}
