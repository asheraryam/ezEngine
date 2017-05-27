﻿
EZ_ALWAYS_INLINE const char* ezWorld::GetName() const
{
  return m_Data.m_sName.GetData();
}

EZ_ALWAYS_INLINE ezUInt32 ezWorld::GetIndex() const
{
  return m_uiIndex;
}

EZ_FORCE_INLINE ezGameObjectHandle ezWorld::CreateObject(const ezGameObjectDesc& desc)
{
  ezGameObject* pNewObject;
  return CreateObject(desc, pNewObject);
}

EZ_FORCE_INLINE bool ezWorld::IsValidObject(const ezGameObjectHandle& object) const
{
  CheckForReadAccess();
  EZ_ASSERT_DEV(object.IsInvalidated() || object.m_InternalId.m_WorldIndex == m_uiIndex,
                "Object does not belong to this world. Expected world id {0} got id {1}", m_uiIndex, object.m_InternalId.m_WorldIndex);

  return m_Data.m_Objects.Contains(object);
}

EZ_FORCE_INLINE bool ezWorld::TryGetObject(const ezGameObjectHandle& object, ezGameObject*& out_pObject)
{
  CheckForReadAccess();
  EZ_ASSERT_DEV(object.IsInvalidated() || object.m_InternalId.m_WorldIndex == m_uiIndex,
                "Object does not belong to this world. Expected world id {0} got id {1}", m_uiIndex, object.m_InternalId.m_WorldIndex);

  return m_Data.m_Objects.TryGetValue(object, out_pObject);
}

EZ_FORCE_INLINE bool ezWorld::TryGetObject(const ezGameObjectHandle& object, const ezGameObject*& out_pObject) const
{
  CheckForReadAccess();
  EZ_ASSERT_DEV(object.IsInvalidated() || object.m_InternalId.m_WorldIndex == m_uiIndex,
    "Object does not belong to this world. Expected world id {0} got id {1}", m_uiIndex, object.m_InternalId.m_WorldIndex);

  ezGameObject* pObject = nullptr;
  bool bResult = m_Data.m_Objects.TryGetValue(object, pObject);
  out_pObject = pObject;
  return bResult;
}

EZ_FORCE_INLINE bool ezWorld::TryGetObjectWithGlobalKey(const ezTempHashedString& sGlobalKey, ezGameObject*& out_pObject)
{
  CheckForReadAccess();
  ezGameObjectId id;
  if (m_Data.m_GlobalKeyToIdTable.TryGetValue(sGlobalKey.GetHash(), id))
  {
    out_pObject = m_Data.m_Objects[id];
    return true;
  }

  return false;
}

EZ_FORCE_INLINE bool ezWorld::TryGetObjectWithGlobalKey(const ezTempHashedString& sGlobalKey, const ezGameObject*& out_pObject) const
{
  CheckForReadAccess();
  ezGameObjectId id;
  if (m_Data.m_GlobalKeyToIdTable.TryGetValue(sGlobalKey.GetHash(), id))
  {
    out_pObject = m_Data.m_Objects[id];
    return true;
  }

  return false;
}

EZ_FORCE_INLINE ezUInt32 ezWorld::GetObjectCount() const
{
  CheckForReadAccess();
  return m_Data.m_ObjectStorage.GetCount();
}

EZ_FORCE_INLINE ezInternal::WorldData::ObjectStorage::Iterator ezWorld::GetObjects()
{
  CheckForWriteAccess();
  return m_Data.m_ObjectStorage.GetIterator(0);
}

EZ_FORCE_INLINE ezInternal::WorldData::ObjectStorage::ConstIterator ezWorld::GetObjects() const
{
  CheckForReadAccess();
  return m_Data.m_ObjectStorage.GetIterator(0);
}

EZ_FORCE_INLINE void ezWorld::Traverse(VisitorFunc visitorFunc, TraversalMethod method /*= DepthFirst*/)
{
  CheckForWriteAccess();

  if (method == DepthFirst)
  {
    m_Data.TraverseDepthFirst(visitorFunc);
  }
  else // method == BreadthFirst
  {
    m_Data.TraverseBreadthFirst(visitorFunc);
  }
}

template <typename ModuleType>
ModuleType* ezWorld::GetOrCreateModule()
{
  CheckForWriteAccess();
  EZ_CHECK_AT_COMPILETIME_MSG(EZ_IS_DERIVED_FROM_STATIC(ezWorldModule, ModuleType),
    "Not a valid module type");

  const ezUInt16 uiTypeId = ModuleType::TypeId();
  if (uiTypeId >= m_Data.m_Modules.GetCount())
  {
    m_Data.m_Modules.SetCount(uiTypeId + 1);
  }

  ModuleType* pModule = static_cast<ModuleType*>(m_Data.m_Modules[uiTypeId]);
  if (pModule == nullptr)
  {
    pModule = EZ_NEW(&m_Data.m_Allocator, ModuleType, this);
    pModule->Initialize();

    m_Data.m_Modules[uiTypeId] = pModule;
    m_Data.m_ModulesToStartSimulation.PushBack(pModule);
  }

  return pModule;
}

template <typename ModuleType>
void ezWorld::DeleteModule()
{
  CheckForWriteAccess();

  const ezUInt16 uiTypeId = ModuleType::TypeId();
  if (uiTypeId < m_Data.m_Modules.GetCount())
  {
    if (ModuleType* pModule = static_cast<ModuleType*>(m_Data.m_Modules[uiTypeId]))
    {
      m_Data.m_Modules[uiTypeId] = nullptr;

      pModule->Deinitialize();
      DeregisterUpdateFunctions(pModule);
      EZ_DELETE(&m_Data.m_Allocator, pModule);
    }
  }
}

template <typename ModuleType>
EZ_FORCE_INLINE ModuleType* ezWorld::GetModule()
{
  CheckForWriteAccess();
  EZ_CHECK_AT_COMPILETIME_MSG(EZ_IS_DERIVED_FROM_STATIC(ezWorldModule, ModuleType),
    "Not a valid module type");

  const ezUInt16 uiTypeId = ModuleType::TypeId();
  if (uiTypeId < m_Data.m_Modules.GetCount())
  {
    return static_cast<ModuleType*>(m_Data.m_Modules[uiTypeId]);
  }

  return nullptr;
}

template <typename ModuleType>
EZ_FORCE_INLINE const ModuleType* ezWorld::GetModule() const
{
  CheckForReadAccess();
  EZ_CHECK_AT_COMPILETIME_MSG(EZ_IS_DERIVED_FROM_STATIC(ezWorldModule, ModuleType),
    "Not a valid module type");

  const ezUInt16 uiTypeId = ModuleType::TypeId();
  if (uiTypeId < m_Data.m_Modules.GetCount())
  {
    return static_cast<const ModuleType*>(m_Data.m_Modules[uiTypeId]);
  }

  return nullptr;
}

template <typename ModuleType>
ModuleType* ezWorld::GetModuleOfBaseType()
{
  CheckForReadAccess();
  EZ_CHECK_AT_COMPILETIME_MSG(EZ_IS_DERIVED_FROM_STATIC(ezWorldModule, ModuleType),
    "Not a valid module type");

  const ezRTTI* pRtti = ezGetStaticRTTI<ModuleType>();
  for (ezWorldModule* pModule : m_Data.m_Modules)
  {
    if (pModule != nullptr && pModule->IsInstanceOf(pRtti))
    {
      return static_cast<ModuleType*>(pModule);
    }
  }

  return nullptr;
}

template <typename ManagerType>
EZ_FORCE_INLINE ManagerType* ezWorld::GetOrCreateComponentManager()
{
  EZ_CHECK_AT_COMPILETIME_MSG(EZ_IS_DERIVED_FROM_STATIC(ezComponentManagerBase, ManagerType),
    "Not a valid component manager type");

  return GetOrCreateModule<ManagerType>();
}

template <typename ManagerType>
EZ_FORCE_INLINE void ezWorld::DeleteComponentManager()
{
  DeleteModule<ManagerType>();
}

template <typename ManagerType>
EZ_FORCE_INLINE ManagerType* ezWorld::GetComponentManager()
{
  EZ_CHECK_AT_COMPILETIME_MSG(EZ_IS_DERIVED_FROM_STATIC(ezComponentManagerBase, ManagerType),
    "Not a valid component manager type");

  return GetModule<ManagerType>();
}

template <typename ManagerType>
EZ_FORCE_INLINE const ManagerType* ezWorld::GetComponentManager() const
{
  EZ_CHECK_AT_COMPILETIME_MSG(EZ_IS_DERIVED_FROM_STATIC(ezComponentManagerBase, ManagerType),
    "Not a valid component manager type");

  return GetModule<ManagerType>();
}

inline bool ezWorld::IsValidComponent(const ezComponentHandle& component) const
{
  CheckForReadAccess();
  const ezUInt16 uiTypeId = component.m_InternalId.m_TypeId;

  if (uiTypeId < m_Data.m_Modules.GetCount())
  {
    if (const ezWorldModule* pModule = m_Data.m_Modules[uiTypeId])
    {
      return static_cast<const ezComponentManagerBase*>(pModule)->IsValidComponent(component);
    }
  }

  return false;
}

template <typename ComponentType>
inline bool ezWorld::TryGetComponent(const ezComponentHandle& component, ComponentType*& out_pComponent)
{
  CheckForWriteAccess();
  EZ_CHECK_AT_COMPILETIME_MSG(EZ_IS_DERIVED_FROM_STATIC(ezComponent, ComponentType), "Not a valid component type");

  const ezUInt16 uiTypeId = component.m_InternalId.m_TypeId;

  if (uiTypeId < m_Data.m_Modules.GetCount())
  {
    if (ezWorldModule* pModule = m_Data.m_Modules[uiTypeId])
    {
      ezComponent* pComponent = nullptr;
      bool bResult = static_cast<ezComponentManagerBase*>(pModule)->TryGetComponent(component, pComponent);
      out_pComponent = ezDynamicCast<ComponentType*>(pComponent);
      return bResult && out_pComponent != nullptr;
    }
  }

  return false;
}

template <typename ComponentType>
inline bool ezWorld::TryGetComponent(const ezComponentHandle& component, const ComponentType*& out_pComponent) const
{
  CheckForReadAccess();
  EZ_CHECK_AT_COMPILETIME_MSG(EZ_IS_DERIVED_FROM_STATIC(ezComponent, ComponentType), "Not a valid component type");

  const ezUInt16 uiTypeId = component.m_InternalId.m_TypeId;

  if (uiTypeId < m_Data.m_Modules.GetCount())
  {
    if (const ezWorldModule* pModule = m_Data.m_Modules[uiTypeId])
    {
      const ezComponent* pComponent = nullptr;
      bool bResult = static_cast<const ezComponentManagerBase*>(pModule)->TryGetComponent(component, pComponent);
      out_pComponent = ezDynamicCast<const ComponentType*>(pComponent);
      return bResult && out_pComponent != nullptr;
    }
  }

  return false;
}

EZ_FORCE_INLINE void ezWorld::SendMessage(const ezGameObjectHandle& receiverObject, ezMessage& msg)
{
  CheckForWriteAccess();

  ezGameObject* pReceiverObject = nullptr;
  if (TryGetObject(receiverObject, pReceiverObject))
  {
    pReceiverObject->SendMessage(msg);
  }
  else
  {
#if EZ_ENABLED(EZ_COMPILE_FOR_DEBUG)
    if (msg.GetDebugMessageRouting())
    {
      ezLog::Warning("ezWorld::SendMessage: The receiver ezGameObject for message of type {0} does not exist.", msg.GetId());
    }
#endif
  }
}

EZ_FORCE_INLINE void ezWorld::SendMessage(const ezComponentHandle& receiverComponent, ezMessage& msg)
{
  CheckForWriteAccess();

  ezComponent* pReceiverComponent = nullptr;
  if (TryGetComponent(receiverComponent, pReceiverComponent))
  {
    pReceiverComponent->SendMessage(msg);
  }
  else
  {
#if EZ_ENABLED(EZ_COMPILE_FOR_DEBUG)
    if (msg.GetDebugMessageRouting())
    {
      ezLog::Warning("ezWorld::SendMessage: The receiver ezComponent for message of type {0} does not exist.", msg.GetId());
    }
#endif
  }
}

EZ_ALWAYS_INLINE void ezWorld::SetWorldSimulationEnabled(bool bEnable)
{
  m_Data.m_bSimulateWorld = bEnable;
}

EZ_ALWAYS_INLINE bool ezWorld::GetWorldSimulationEnabled() const
{
  return m_Data.m_bSimulateWorld;
}

EZ_ALWAYS_INLINE ezTask* ezWorld::GetUpdateTask()
{
  return &m_UpdateTask;
}

EZ_FORCE_INLINE ezSpatialSystem& ezWorld::GetSpatialSystem()
{
  CheckForWriteAccess();

  return *(m_Data.m_pSpatialSystem.Borrow());
}

EZ_FORCE_INLINE const ezSpatialSystem& ezWorld::GetSpatialSystem() const
{
  CheckForReadAccess();

  return *(m_Data.m_pSpatialSystem.Borrow());
}

EZ_ALWAYS_INLINE void ezWorld::GetCoordinateSystem(const ezVec3& vGlobalPosition, ezCoordinateSystem& out_CoordinateSystem) const
{
  m_Data.m_pCoordinateSystemProvider->GetCoordinateSystem(vGlobalPosition, out_CoordinateSystem);
}

EZ_FORCE_INLINE void ezWorld::SetCoordinateSystemProvider(ezUniquePtr<ezCoordinateSystemProvider>&& pProvider)
{
  EZ_ASSERT_DEV(pProvider != nullptr, "Coordinate System Provider must not be null");

  m_Data.m_pCoordinateSystemProvider = std::move(pProvider);
  m_Data.m_pCoordinateSystemProvider->m_pOwnerWorld = this;
}

EZ_ALWAYS_INLINE ezCoordinateSystemProvider& ezWorld::GetCoordinateSystemProvider()
{
  return *(m_Data.m_pCoordinateSystemProvider.Borrow());
}

EZ_ALWAYS_INLINE const ezCoordinateSystemProvider& ezWorld::GetCoordinateSystemProvider() const
{
  return *(m_Data.m_pCoordinateSystemProvider.Borrow());
}

EZ_ALWAYS_INLINE ezClock& ezWorld::GetClock()
{
  return m_Data.m_Clock;
}

EZ_ALWAYS_INLINE const ezClock& ezWorld::GetClock() const
{
  return m_Data.m_Clock;
}

EZ_ALWAYS_INLINE ezRandom& ezWorld::GetRandomNumberGenerator()
{
  return m_Data.m_Random;
}

EZ_ALWAYS_INLINE ezAllocatorBase* ezWorld::GetAllocator()
{
  return &m_Data.m_Allocator;
}

EZ_ALWAYS_INLINE ezInternal::WorldLargeBlockAllocator* ezWorld::GetBlockAllocator()
{
  return &m_Data.m_BlockAllocator;
}

EZ_ALWAYS_INLINE ezInternal::WorldData::ReadMarker& ezWorld::GetReadMarker() const
{
  return m_Data.m_ReadMarker;
}

EZ_ALWAYS_INLINE ezInternal::WorldData::WriteMarker& ezWorld::GetWriteMarker()
{
  return m_Data.m_WriteMarker;
}

EZ_FORCE_INLINE void ezWorld::SetUserData(void* pUserData)
{
  CheckForWriteAccess();

  m_Data.m_pUserData = pUserData;
}

EZ_FORCE_INLINE void* ezWorld::GetUserData() const
{
  CheckForReadAccess();

  return m_Data.m_pUserData;
}

//static
EZ_ALWAYS_INLINE ezUInt32 ezWorld::GetWorldCount()
{
  return s_Worlds.GetCount();
}

//static
EZ_ALWAYS_INLINE ezWorld* ezWorld::GetWorld(ezUInt32 uiIndex)
{
  return s_Worlds[uiIndex];
}

EZ_FORCE_INLINE void ezWorld::CheckForReadAccess() const
{
  EZ_ASSERT_DEV(m_Data.m_iReadCounter > 0, "Trying to read from World '{0}', but it is not marked for reading.", GetName());
}

EZ_FORCE_INLINE void ezWorld::CheckForWriteAccess() const
{
  EZ_ASSERT_DEV(m_Data.m_WriteThreadID == ezThreadUtils::GetCurrentThreadID(), "Trying to write to World '{0}', but it is not marked for writing.", GetName());
}

EZ_ALWAYS_INLINE ezGameObject* ezWorld::GetObjectUnchecked(ezUInt32 uiIndex) const
{
  return m_Data.m_Objects.GetValueUnchecked(uiIndex);
}
