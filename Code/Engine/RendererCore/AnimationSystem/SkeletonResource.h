#pragma once

#include <Core/ResourceManager/Resource.h>
#include <RendererCore/AnimationSystem/Declarations.h>
#include <RendererCore/AnimationSystem/Skeleton.h>
#include <RendererCore/RendererCoreDLL.h>

struct ezSkeletonResourceGeometry
{
  // scale is used to resize a unit sphere / box / capsule
  ezTransform m_Transform;
  ezUInt16 m_uiAttachedToJoint = 0;
  ezEnum<ezSkeletonJointGeometryType> m_Type;
};

struct EZ_RENDERERCORE_DLL ezSkeletonResourceDescriptor
{
  ezSkeletonResourceDescriptor();
  ~ezSkeletonResourceDescriptor();
  ezSkeletonResourceDescriptor(const ezSkeletonResourceDescriptor& rhs) = delete;
  ezSkeletonResourceDescriptor(ezSkeletonResourceDescriptor&& rhs);
  void operator=(ezSkeletonResourceDescriptor&& rhs);
  void operator=(const ezSkeletonResourceDescriptor& rhs) = delete;

  ezResult Serialize(ezStreamWriter& stream) const;
  ezResult Deserialize(ezStreamReader& stream);

  ezSkeleton m_Skeleton;
  ezDynamicArray<ezSkeletonResourceGeometry> m_Geometry;
};

typedef ezTypedResourceHandle<class ezSkeletonResource> ezSkeletonResourceHandle;

class EZ_RENDERERCORE_DLL ezSkeletonResource : public ezResource
{
  EZ_ADD_DYNAMIC_REFLECTION(ezSkeletonResource, ezResource);
  EZ_RESOURCE_DECLARE_COMMON_CODE(ezSkeletonResource);
  EZ_RESOURCE_DECLARE_CREATEABLE(ezSkeletonResource, ezSkeletonResourceDescriptor);

public:
  ezSkeletonResource();
  ~ezSkeletonResource();

  const ezSkeletonResourceDescriptor& GetDescriptor() const { return m_Descriptor; }

private:
  virtual ezResourceLoadDesc UnloadData(Unload WhatToUnload) override;
  virtual ezResourceLoadDesc UpdateContent(ezStreamReader* Stream) override;
  virtual void UpdateMemoryUsage(MemoryUsage& out_NewMemoryUsage) override;

  ezSkeletonResourceDescriptor m_Descriptor;
};
