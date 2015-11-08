#pragma once

#include <ToolsFoundation/Basics.h>
#include <ToolsFoundation/Command/Command.h>
#include <ToolsFoundation/Document/Document.h>

class ezDuplicateObjectsCommand : public ezCommand
{
  EZ_ADD_DYNAMIC_REFLECTION(ezDuplicateObjectsCommand, ezCommand);

public:
  ezDuplicateObjectsCommand();

public: // Properties
  ezString m_sJsonGraph;
  ezString m_sParentNodes;

private:
  virtual ezStatus Do(bool bRedo) override;
  virtual ezStatus Undo(bool bFireEvents) override;
  virtual void Cleanup(CommandState state) override;

private:
  struct DuplicatedObject
  {
    ezDocumentObject* m_pObject;
    ezDocumentObject* m_pParent;
    ezString m_sParentProperty;
    ezVariant m_Index;
  };

  ezHybridArray<DuplicatedObject, 4> m_DuplicatedObjects;
};



