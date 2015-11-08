#pragma once

#include <EditorFramework/DocumentWindow3D/EditorInputContext.h>
#include <EditorFramework/Gizmos/GizmoBase.h>
#include <QPoint>

class QWidget;
class ezCamera;

class ezOrthoGizmoContext : public ezEditorInputContext
{
  EZ_ADD_DYNAMIC_REFLECTION(ezOrthoGizmoContext, ezEditorInputContext);

public:

  ezOrthoGizmoContext(ezQtEngineDocumentWindow* pOwnerWindow, ezQtEngineViewWidget* pOwnerView, const ezCamera* pCamera);

  void SetWindowConfig(const ezVec2I32& viewport)
  {
    m_Viewport = viewport;
  }

  virtual void FocusLost(bool bCancel);

  virtual ezEditorInut mousePressEvent(QMouseEvent* e) override;
  virtual ezEditorInut mouseReleaseEvent(QMouseEvent* e) override;
  virtual ezEditorInut mouseMoveEvent(QMouseEvent* e) override;

  ezEvent<const ezGizmo::GizmoEvent&> m_GizmoEvents;

  const ezVec3& GetTranslationResult() const { return m_vTranslationResult; }
  const ezQuat& GetRotationResult() const { return m_qRotationResult; }
  float GetScalingResult() const { return m_fScalingResult; }

protected:
  virtual void OnSetOwner(ezQtEngineDocumentWindow* pOwnerWindow, ezQtEngineViewWidget* pOwnerView) override {}

private:
  bool IsViewInOthoMode() const;

  QPoint m_LastMousePos;
  ezVec3 m_vTranslationResult;
  ezQuat m_qRotationResult;
  float m_fScaleMouseMove;
  float m_fScalingResult;
  bool m_bCanInteract;
  const ezCamera* m_pCamera;
  ezVec2I32 m_Viewport;
};

