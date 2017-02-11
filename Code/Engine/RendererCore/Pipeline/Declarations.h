#pragma once

#include <Foundation/Reflection/Reflection.h>
#include <Core/World/Declarations.h>
#include <RendererCore/Basics.h>
#include <RendererFoundation/Basics.h>
#include <RendererFoundation/Descriptors/Descriptors.h>

class ezCamera;
class ezExtractedRenderData;
class ezExtractor;
class ezView;
class ezRenderData;
class ezRenderDataBatch;
class ezRenderPipeline;
class ezRenderPipelinePass;
class ezGALContext;
class ezRenderContext;
class ezDebugRendererContext;

struct ezNodePin;
struct ezRenderPipelinePassConnection;
struct ezViewData;

struct ezRenderViewContext
{
  const ezCamera* m_pCamera;
  const ezViewData* m_pViewData;
  ezRenderContext* m_pRenderContext;
  
  const ezDebugRendererContext* m_pWorldDebugContext;
  const ezDebugRendererContext* m_pViewDebugContext;
};

/// \brief This is the base class for types that handle rendering of different object types.
///
/// E.g. there are different renderers for meshes, particle effects, light sources, etc.
class EZ_RENDERERCORE_DLL ezRenderer : public ezReflectedClass
{
  EZ_ADD_DYNAMIC_REFLECTION(ezRenderer, ezReflectedClass);

public:
  virtual void GetSupportedRenderDataTypes(ezHybridArray<const ezRTTI*, 8>& types) = 0;

  virtual void RenderBatch(const ezRenderViewContext& renderViewContext, ezRenderPipelinePass* pPass, const ezRenderDataBatch& batch) = 0;
};

/// \brief Usage hint of a camera/view.
struct EZ_RENDERERCORE_DLL ezCameraUsageHint
{
  typedef ezInt8 StorageType;

  enum Enum
  {
    None,
    MainView,
    EditorView,
    Thumbnail,

    ENUM_COUNT,

    Default = None,
  };
};

EZ_DECLARE_REFLECTABLE_TYPE(EZ_RENDERERCORE_DLL, ezCameraUsageHint);
