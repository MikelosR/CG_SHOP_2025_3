#ifndef EXTRA_GRAPHICS_H
#define EXTRA_GRAPHICS_H

#include <QtWidgets/QApplication>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QGraphicsScene>
#include <QtGui/QPainter>
#include <sstream>
#include <QtWidgets/QGraphicsTextItem>
#include <QScrollBar>
#include <QtWidgets/QToolTip>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include "functions.h"

class CDTGraphicsView : public QGraphicsView {
public:
    CDTGraphicsView(CDT& cdt, Polygon& polygon, QWidget* parent = nullptr) 
        : QGraphicsView(parent), cdt(cdt), polygon(polygon) {
        //Initialize the scene and set it for this view
        QGraphicsScene* scene = new QGraphicsScene(this);
        setScene(scene);
        //Draw the triangulation
        drawTriangulation();
        //Fit the entire triangulation into the view window
        fitInView(scene->sceneRect(), Qt::KeepAspectRatio);
    }    

protected:
    void drawTriangulation() {
        QGraphicsScene* scene = this->scene();

        // Calculate the bounding box for the points to set the scene rectangle
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double max_y = std::numeric_limits<double>::lowest();

        for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
            Point_2 p = vit->point();
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            min_x = std::min(min_x, x);
            min_y = std::min(min_y, y);
            max_x = std::max(max_x, x);
            max_y = std::max(max_y, y);
        }

        // Include the polygon bounding box in the scene bounds
        double polygon_min_x = std::numeric_limits<double>::max();
        double polygon_min_y = std::numeric_limits<double>::max();
        double polygon_max_x = std::numeric_limits<double>::lowest();
        double polygon_max_y = std::numeric_limits<double>::lowest();
        for (auto vertex = polygon.vertices_begin(); vertex != polygon.vertices_end(); ++vertex) {
            double x = CGAL::to_double(vertex->x());
            double y = CGAL::to_double(vertex->y());
            polygon_min_x = std::min(polygon_min_x, x);
            polygon_min_y = std::min(polygon_min_y, y);
            polygon_max_x = std::max(polygon_max_x, x);
            polygon_max_y = std::max(polygon_max_y, y);
        }

        double margin = 30; // Padding around the triangulation
        double sceneHeight = max_y - min_y;
        scene->setSceneRect(std::min(min_x, polygon_min_x) - margin, 
                            std::min(min_y, polygon_min_y) - margin, 
                            std::max(max_x, polygon_max_x) - min_x + 2 * margin, 
                            std::max(max_y, polygon_max_y) - min_y + 2 * margin);

        // --- Draw the Polygon Boundary ---
        QPen polygonPen(Qt::blue); // Choose a color for the polygon boundary
        QPolygonF polygonPoints;
        for (auto vertex = polygon.vertices_begin(); vertex != polygon.vertices_end(); ++vertex) {
            double x = CGAL::to_double(vertex->x());
            double y = sceneHeight - CGAL::to_double(vertex->y()) + min_y;
            polygonPoints << QPointF(x, y);
        }
        scene->addPolygon(polygonPoints, polygonPen);

        // --- Draw the Edges of the Triangulation ---
        QPen edgePen(QColor(0, 0, 0, 255)); // Fully opaque black color
        edgePen.setWidth(1); // Optional: thicker edges for visibility

        for (auto eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit) {
            auto source = eit->first->vertex((eit->second + 1) % 3)->point();
            auto target = eit->first->vertex((eit->second + 2) % 3)->point();
            scene->addLine(
                CGAL::to_double(source.x()), sceneHeight - CGAL::to_double(source.y()) + min_y,
                CGAL::to_double(target.x()), sceneHeight - CGAL::to_double(target.y()) + min_y, edgePen);
        }

        // --- Draw the Vertices of the Triangulation ---
        QBrush regularVertexBrush(QColor(255, 0, 0, 255)); // Fully opaque red color
        QPen vertexPen(Qt::black); // Optional: black border for better visibility       

        for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
            Point_2 p = vit->point();
            double x = CGAL::to_double(p.x());
            double y = sceneHeight - CGAL::to_double(p.y()) + min_y;

            // Draw the point red color
            scene->addEllipse(x - 10, y - 10, 20, 20, vertexPen, regularVertexBrush);

            // Add coordinates as text
            std::ostringstream oss;
            oss << "(" << x << ", " << CGAL::to_double(p.y()) << ")";
            auto textItem = scene->addText(QString::fromStdString(oss.str()));
            textItem->setPos(x + 15, y + 15);

            textItem->setToolTip(QString::fromStdString(oss.str()));

            // Set the font size
            QFont font = textItem->font();
            font.setPointSize(25);
            textItem->setFont(font);
        }

        // --- Draw Obtuse Triangles and Vertices ---
        QBrush obtuseTriangleBrush(QColor(255, 0, 0, 120)); // Darker red for obtuse triangles
        QBrush obtuseVertexBrush(QColor(0, 255, 0, 255));   // Green for obtuse vertices inside the boundary
        QBrush obtuseOutsideBrush(QColor(0, 0, 0, 128)); // Light black for obtuse triangles outside polygon

        for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            Point_2 p1 = fit->vertex(0)->point();
            Point_2 p2 = fit->vertex(1)->point();
            Point_2 p3 = fit->vertex(2)->point();

            // Check if triangle is obtuse and inside the region boundary
            bool isObtuseTriangle = is_obtuse(fit);
            bool faceInside = is_face_inside_region(fit, polygon);

            // Check if the triangle is obtuse
            if (isObtuseTriangle && faceInside) {
                // Draw the triangle with a light color
                QPolygonF triangle;
                triangle << QPointF(CGAL::to_double(p1.x()), sceneHeight - CGAL::to_double(p1.y()) + min_y)
                        << QPointF(CGAL::to_double(p2.x()), sceneHeight - CGAL::to_double(p2.y()) + min_y)
                        << QPointF(CGAL::to_double(p3.x()), sceneHeight - CGAL::to_double(p3.y()) + min_y);
                scene->addPolygon(triangle, QPen(Qt::NoPen), obtuseTriangleBrush);

                // Identify and draw vertices with the correct color
                for (int i = 0; i < 3; ++i) {
                    Point_2 vertex = fit->vertex(i)->point();
                    bool isObtuseVertex = (find_obtuse_vertex(p1, p2, p3) == vertex);
                    double x = CGAL::to_double(vertex.x());
                    double y = sceneHeight - CGAL::to_double(vertex.y()) + min_y;

                    // Draw vertices based on their color category
                    if (polygon.bounded_side(vertex) != CGAL::ON_UNBOUNDED_SIDE) {
                        // If the vertex is obtuse, draw it green
                        if (isObtuseVertex) {
                            QBrush vertexBrush = obtuseVertexBrush;
                            scene->addEllipse(x - 10, y - 10, 20, 20, QPen(Qt::NoPen), vertexBrush);
                        }
                    }
                }
            }

            // If face is outside of the boundary
            if (!faceInside) {
                // Draw the triangle with a light black color
                QPolygonF triangle2;
                triangle2 << QPointF(CGAL::to_double(p1.x()), sceneHeight - CGAL::to_double(p1.y()) + min_y)
                        << QPointF(CGAL::to_double(p2.x()), sceneHeight - CGAL::to_double(p2.y()) + min_y)
                        << QPointF(CGAL::to_double(p3.x()), sceneHeight - CGAL::to_double(p3.y()) + min_y);
                scene->addPolygon(triangle2, QPen(Qt::NoPen), obtuseOutsideBrush);
            }
        }
    }


    //Override mouse press event to start the dragging
    void mousePressEvent(QMouseEvent* event) override {
        if (event->button() == Qt::LeftButton) {
            dragStartPos = event->pos(); //Store the initial position
            dragInProgress = true; //Start dragging
        }
        QGraphicsView::mousePressEvent(event); //Call base class handler
    }

    //Override mouse move event to handle the dragging and tooltip display
    void mouseMoveEvent(QMouseEvent* event) override {
        //Handle dragging
        if (dragInProgress) {
            QPointF delta = event->pos() - dragStartPos; //Calculate the movement
            horizontalScrollBar()->setValue(horizontalScrollBar()->value() - delta.x());
            verticalScrollBar()->setValue(verticalScrollBar()->value() - delta.y());
            dragStartPos = event->pos(); //Update the position
        } else {
            //Handle tooltip display
            QPointF scenePos = mapToScene(event->pos());
            QGraphicsItem* item = scene()->itemAt(scenePos, QTransform());

            if (item) {
                //Assuming you label your points with their coordinates or an identifier
                QString tooltipText = item->toolTip(); //Get the tooltip
                QToolTip::showText(event->globalPos(), tooltipText); //Show the tooltip at the cursor position
            } else {
                QToolTip::hideText(); //Hide the tooltip if not over an item
            }
        }
        QGraphicsView::mouseMoveEvent(event); //Call base class handler
    }


    // Override mouse release event to end the dragging
    void mouseReleaseEvent(QMouseEvent* event) override {
        if (event->button() == Qt::LeftButton) {
            dragInProgress = false; //End dragging
        }
        QGraphicsView::mouseReleaseEvent(event); //Call base class handler
    }

    // Override wheel event for zooming
    void wheelEvent(QWheelEvent* event) override {
        if (event->angleDelta().y() > 0) {
            scale(1.2, 1.2); // Zoom in
        } else {
            scale(1 / 1.2, 1 / 1.2); // Zoom out
        }
        event->accept();
    }

private:
    CDT& cdt;
    const Polygon& polygon;
    QPointF dragStartPos; //To store the position when dragging starts
    bool dragInProgress = false; //To track if dragging is happening
};

#endif //EXTRA_GRAPHICS_H