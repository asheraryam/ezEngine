#pragma once

#include <Foundation/Types/Id.h>

namespace ezModelImporter
{
  typedef ezGenericId<24, 8> ObjectId;

  /// A handle to a hierarchy object.
  /// \see ezModelImporter::Scene
  class ObjectHandle
  {
  public:
    EZ_DECLARE_POD_TYPE();

    // Cannot use EZ_DECLARE_HANDLE_TYPE since type is an integral part.

    ObjectHandle() : m_Type(NULLREF), m_Id() {}

    void operator = (const ObjectHandle& rhs);
    bool operator == (const ObjectHandle& rhs) const;

    bool IsValid() const { return m_Type != NULLREF && m_Id != ObjectId(); }

    enum Type
    {
      NODE,
      MESH,
      //LIGHT,
      //CAMERA,

      NULLREF = -1  // null reference.
    };

    Type GetType() const { return m_Type; }

  private:
    // Only the scene is allowed to create these handles.
    friend class Scene;
    // Access for easier hashing.
    friend struct ezHashHelper<ObjectHandle>;


    ObjectHandle(Type type, ObjectId id) : m_Type(type), m_Id(id) {}


    Type m_Type;
    ObjectId m_Id;
  };

  /// Hash function for ObjectHandle
  template <>
  struct ezHashHelper<ObjectHandle>
  {
    static ezUInt32 Hash(ObjectHandle value)
    {
      // First 24 bit are instance. It's unlikely those will ever be fully used.
      // Type takes exactly 2 bit as of writing.
      ezUInt32 idValue = value.m_Id.m_Data; 
      return idValue ^ (value.m_Type << 22);
    }

    static bool Equal(ObjectHandle a, ObjectHandle b)
    {
      return a == b;
    }
  };

  typedef ezGenericId<22, 10> MaterialId;

  /// A handle to a material.
  /// \see ezModelImporter::Scene
  class MaterialHandle
  {
    EZ_DECLARE_HANDLE_TYPE(MaterialHandle, MaterialId);
  };

}
