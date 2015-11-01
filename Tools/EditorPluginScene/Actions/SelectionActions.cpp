#include <PCH.h>
#include <GuiFoundation/Action/ActionManager.h>
#include <GuiFoundation/Action/ActionMapManager.h>
#include <EditorPluginScene/Actions/SelectionActions.h>
#include <EditorPluginScene/Scene/SceneDocument.h>
#include <QFileDialog>

EZ_BEGIN_DYNAMIC_REFLECTED_TYPE(ezSelectionAction, ezButtonAction, 0, ezRTTINoAllocator);
EZ_END_DYNAMIC_REFLECTED_TYPE();

ezActionDescriptorHandle ezSelectionActions::s_hSelectionCategory;
ezActionDescriptorHandle ezSelectionActions::s_hShowInScenegraph;
ezActionDescriptorHandle ezSelectionActions::s_hFocusOnSelection;
ezActionDescriptorHandle ezSelectionActions::s_hFocusOnSelectionAllViews;
ezActionDescriptorHandle ezSelectionActions::s_hGroupSelectedItems;
ezActionDescriptorHandle ezSelectionActions::s_hHideSelectedObjects;
ezActionDescriptorHandle ezSelectionActions::s_hHideUnselectedObjects;
ezActionDescriptorHandle ezSelectionActions::s_hShowHiddenObjects;
ezActionDescriptorHandle ezSelectionActions::s_hCreatePrefab;

void ezSelectionActions::RegisterActions()
{
  s_hSelectionCategory = EZ_REGISTER_CATEGORY("SelectionCategory");
  s_hShowInScenegraph = EZ_REGISTER_ACTION_1("ActionShowInScenegraph", ezActionScope::Document, "Scene - Selection", "Ctrl+T", ezSelectionAction, ezSelectionAction::ActionType::ShowInScenegraph);
  s_hFocusOnSelection = EZ_REGISTER_ACTION_1("ActionFocusOnSelection", ezActionScope::Document, "Scene - Selection", "F", ezSelectionAction, ezSelectionAction::ActionType::FocusOnSelection);
  s_hFocusOnSelectionAllViews = EZ_REGISTER_ACTION_1("ActionFocusOnSelectionAllViews", ezActionScope::Document, "Scene - Selection", "Ctrl+K,Ctrl+F", ezSelectionAction, ezSelectionAction::ActionType::FocusOnSelectionAllViews);
  s_hGroupSelectedItems = EZ_REGISTER_ACTION_1("ActionGroupSelectedItems", ezActionScope::Document, "Scene - Selection", "G", ezSelectionAction, ezSelectionAction::ActionType::GroupSelectedItems);
  s_hHideSelectedObjects = EZ_REGISTER_ACTION_1("ActionHideSelectedObjects", ezActionScope::Document, "Scene - Selection", "H", ezSelectionAction, ezSelectionAction::ActionType::HideSelectedObjects);
  s_hHideUnselectedObjects = EZ_REGISTER_ACTION_1("ActionHideUnselectedObjects", ezActionScope::Document, "Scene - Selection", "Shift+H", ezSelectionAction, ezSelectionAction::ActionType::HideUnselectedObjects);
  s_hShowHiddenObjects = EZ_REGISTER_ACTION_1("ActionShowHiddenObjects", ezActionScope::Document, "Scene - Selection", "Ctrl+H", ezSelectionAction, ezSelectionAction::ActionType::ShowHiddenObjects);
  s_hCreatePrefab = EZ_REGISTER_ACTION_1("ActionCreatePrefab", ezActionScope::Document, "Scene - Selection", "Ctrl+P, Ctrl+C", ezSelectionAction, ezSelectionAction::ActionType::CreatePrefab);
}

void ezSelectionActions::UnregisterActions()
{
  ezActionManager::UnregisterAction(s_hSelectionCategory);
  ezActionManager::UnregisterAction(s_hShowInScenegraph);
  ezActionManager::UnregisterAction(s_hFocusOnSelection);
  ezActionManager::UnregisterAction(s_hFocusOnSelectionAllViews);
  ezActionManager::UnregisterAction(s_hGroupSelectedItems);
  ezActionManager::UnregisterAction(s_hHideSelectedObjects);
  ezActionManager::UnregisterAction(s_hHideUnselectedObjects);
  ezActionManager::UnregisterAction(s_hShowHiddenObjects);
  ezActionManager::UnregisterAction(s_hCreatePrefab);
}

void ezSelectionActions::MapActions(const char* szMapping, const char* szPath)
{
  ezActionMap* pMap = ezActionMapManager::GetActionMap(szMapping);
  EZ_ASSERT_DEV(pMap != nullptr, "The given mapping ('%s') does not exist, mapping the actions failed!", szMapping);

  ezStringBuilder sSubPath(szPath, "/SelectionCategory");

  pMap->MapAction(s_hSelectionCategory, szPath, 5.0f);
  
  pMap->MapAction(s_hShowInScenegraph, sSubPath, 2.0f);
  pMap->MapAction(s_hFocusOnSelection, sSubPath, 3.0f);
  pMap->MapAction(s_hFocusOnSelectionAllViews, sSubPath, 3.5f);
  pMap->MapAction(s_hGroupSelectedItems, sSubPath, 3.7f);
  pMap->MapAction(s_hHideSelectedObjects, sSubPath, 4.0f);
  pMap->MapAction(s_hHideUnselectedObjects, sSubPath, 5.0f);
  pMap->MapAction(s_hShowHiddenObjects, sSubPath, 6.0f);
}


void ezSelectionActions::MapContextMenuActions(const char* szMapping, const char* szPath)
{
  ezActionMap* pMap = ezActionMapManager::GetActionMap(szMapping);
  EZ_ASSERT_DEV(pMap != nullptr, "The given mapping ('%s') does not exist, mapping the actions failed!", szMapping);

  ezStringBuilder sSubPath(szPath, "/SelectionCategory");

  pMap->MapAction(s_hSelectionCategory, szPath, 5.0f);
  
  pMap->MapAction(s_hFocusOnSelectionAllViews, sSubPath, 1.0f);
  pMap->MapAction(s_hGroupSelectedItems, sSubPath, 2.0f);
  pMap->MapAction(s_hHideSelectedObjects, sSubPath, 3.0f);
  pMap->MapAction(s_hCreatePrefab, sSubPath, 4.0f);
}

ezSelectionAction::ezSelectionAction(const ezActionContext& context, const char* szName, ezSelectionAction::ActionType type) : ezButtonAction(context, szName, false, "")
{
  m_Type = type;
  m_pSceneDocument = static_cast<ezSceneDocument*>(context.m_pDocument);

  switch (m_Type)
  {
  case ActionType::ShowInScenegraph:
    SetIconPath(":/GuiFoundation/Icons/Scenegraph16.png");
    break;
  case ActionType::FocusOnSelection:
    SetIconPath(":/GuiFoundation/Icons/FocusOnSelection16.png");
    break;
  case ActionType::FocusOnSelectionAllViews:
    SetIconPath(":/GuiFoundation/Icons/FocusOnSelection16.png"); /// \todo Icon
    break;
  case ActionType::GroupSelectedItems:
    SetIconPath(":/GuiFoundation/Icons/GroupSelection16.png");
    break;
  case ActionType::HideSelectedObjects:
    SetIconPath(":/GuiFoundation/Icons/HideSelected16.png");
    break;
  case ActionType::HideUnselectedObjects:
    SetIconPath(":/GuiFoundation/Icons/HideUnselected16.png");
    break;
  case ActionType::ShowHiddenObjects:
    SetIconPath(":/GuiFoundation/Icons/ShowHidden16.png");
    break;
  case ActionType::CreatePrefab:
    //SetIconPath(":/GuiFoundation/Icons/ShowHidden16.png"); /// \todo icon
    break;
  }

  UpdateEnableState();

  m_Context.m_pDocument->GetSelectionManager()->m_Events.AddEventHandler(ezMakeDelegate(&ezSelectionAction::SelectionEventHandler, this));
}


ezSelectionAction::~ezSelectionAction()
{
  m_Context.m_pDocument->GetSelectionManager()->m_Events.RemoveEventHandler(ezMakeDelegate(&ezSelectionAction::SelectionEventHandler, this));
}

void ezSelectionAction::Execute(const ezVariant& value)
{
  switch (m_Type)
  {
  case ActionType::ShowInScenegraph:
    m_pSceneDocument->TriggerShowSelectionInScenegraph();
    return;
  case ActionType::FocusOnSelection:
    m_pSceneDocument->TriggerFocusOnSelection(false);
    return;
  case ActionType::FocusOnSelectionAllViews:
    m_pSceneDocument->TriggerFocusOnSelection(true);
    return;
  case ActionType::GroupSelectedItems:
    m_pSceneDocument->GroupSelection();
    return;
  case ActionType::HideSelectedObjects:
    m_pSceneDocument->ShowOrHideSelectedObjects(ezSceneDocument::ShowOrHide::Hide);
    break;
  case ActionType::HideUnselectedObjects:
    m_pSceneDocument->HideUnselectedObjects();
    break;
  case ActionType::ShowHiddenObjects:
    m_pSceneDocument->ShowOrHideAllObjects(ezSceneDocument::ShowOrHide::Show);
    break;
  case ActionType::CreatePrefab:
    {
      ezStringBuilder sFile = QFileDialog::getSaveFileName(QApplication::activeWindow(), QLatin1String("Create Prefab"), "", QString::fromUtf8("*.ezPrefab")).toUtf8().data();

      if (!sFile.IsEmpty())
      {
        sFile.ChangeFileExtension("ezPrefab");
        
        auto res = m_pSceneDocument->CreatePrefabDocumentFromSelection(sFile);
        ezUIServices::MessageBoxStatus(res, "Failed to create Prefab", "Successfully created Prefab");
      }
    }
    return;
  }
}

void ezSelectionAction::SelectionEventHandler(const ezSelectionManager::Event& e)
{
  UpdateEnableState();

}

void ezSelectionAction::UpdateEnableState()
{
  if (m_Type == ActionType::FocusOnSelection ||
      m_Type == ActionType::FocusOnSelectionAllViews ||
      m_Type == ActionType::HideSelectedObjects ||
      m_Type == ActionType::ShowInScenegraph ||
      m_Type == ActionType::HideUnselectedObjects)
  {
    SetEnabled(!m_Context.m_pDocument->GetSelectionManager()->IsSelectionEmpty());
  }

  if (m_Type == ActionType::GroupSelectedItems)
  {
    SetEnabled(m_Context.m_pDocument->GetSelectionManager()->GetSelection().GetCount() > 1);
  }

  if (m_Type == ActionType::CreatePrefab)
  {
    SetEnabled(m_Context.m_pDocument->GetSelectionManager()->GetSelection().GetCount() == 1);
  }
}

