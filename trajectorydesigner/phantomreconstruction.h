#ifndef PHANTOMRECONSTRUCTION_H
#define PHANTOMRECONSTRUCTION_H

#include <QWidget>

#include <memory>

QT_FORWARD_DECLARE_CLASS(Phantom)

class PhantomReconstruction : public QWidget
{
	Q_OBJECT
public:
	explicit PhantomReconstruction(QWidget *parent = 0);

signals:

public slots:
private:
	std::shared_ptr<Phantom> m_phantom;
};

#endif // PHANTOMRECONSTRUCTION_H
