
#include "Camera_manip_zoom.h"

#include "Message_manager.h"


Camera_manip_zoom::Camera_manip_zoom(Camera& camera) : 
  Camera_manip(camera)
{
}

void Camera_manip_zoom::mouse_move_event(QMouseEvent* e)
{
  if (m_middle_mouse_button_down)
  {
    const float zoom_scale_factor = 0.01f;
    const auto distance = zoom_scale_factor * m_diff.y();
    m_camera.move_forward(distance);
  }
}
void Camera_manip_zoom::mouse_release_event(QMouseEvent* e)
{
  if (e->button() == Qt::MiddleButton)
  {
    Message_manager::notify_all("zoom_changed");
  }
}

